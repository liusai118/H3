import os
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from torch.amp import autocast
from tqdm.auto import tqdm

# -------------------------
# Repro / speed controls
# -------------------------
def seed_everything(seed=123):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


torch.backends.cudnn.benchmark = True
seed_everything(123)


try:
    torch.set_float32_matmul_precision("high")
except Exception:
    pass

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
use_amp = (device.type == "cuda")
amp_device = "cuda" if use_amp else "cpu"

print("Device:", device)
if device.type == "cuda":
    print("GPU:", torch.cuda.get_device_name(0))

# -------------------------
# Config (EDIT PATHS)
# -------------------------
TEST_CSV = "./pos.csv"
SEQ_COL = "sequence"

WEIGHTS_ROOT = "./out"   # rep1..rep10/best_model.pt
OUT_DIR = "./pos_pre"
os.makedirs(OUT_DIR, exist_ok=True)



REP_THRESH = {1: 0.44, 2: 0.66, 3: 0.55, 4: 0.76, 5: 0.60, 6: 0.32, 7: 0.53, 8: 0.57, 9: 0.46, 10: 0.50}

SEQ_LEN = 256
N_FILL = 0.25

BATCH_SIZE = 512
NUM_WORKERS = 0
PIN_MEMORY = (device.type == "cuda")

IG_BATCH = 64
TOP_PCT_BY_FINAL = 1.0
IG_STEPS = 50
WINDOW = 10
TOPK_WINDOWS = 3
EDGE_TRIM = 5
BASELINE_VALUE = 0.25

RUN_PROB_ONLY = True
RUN_WITH_IG = True


_lut = np.full(256, -1, dtype=np.int16)
for ch, idx in [("A", 0), ("C", 1), ("G", 2), ("T", 3)]:
    _lut[ord(ch)] = idx
    _lut[ord(ch.lower())] = idx

def one_hot_encode_fast_with_len(seq: str, L: int = 256, n_fill: float = 0.25):
    if not isinstance(seq, str):
        seq = str(seq)
    true_len = min(len(seq), L)
    if len(seq) >= L:
        s = seq[:L]
    else:
        s = seq + ("N" * (L - len(seq)))
    b = np.frombuffer(s.encode("ascii", "replace"), dtype=np.uint8)
    if b.size < L:
        b = np.pad(b, (0, L - b.size), constant_values=ord("N"))
    elif b.size > L:
        b = b[:L]
    idx = _lut[b]
    x = np.full((4, L), n_fill, dtype=np.float32)
    pos = np.where(idx >= 0)[0]
    if pos.size > 0:
        x[:, pos] = 0.0
        x[idx[pos], pos] = 1.0
    return x, int(true_len)

def encode_sequences(sequences: list[str], seq_len: int, n_fill: float):
    tmp = [one_hot_encode_fast_with_len(s, seq_len, n_fill=n_fill)
           for s in tqdm(sequences, desc="Encoding TEST", unit="seq")]
    X = np.stack([t[0] for t in tmp]).astype(np.float32)
    lens = np.array([t[1] for t in tmp], dtype=np.int64)
    return X, lens

# -------------------------
# -------------------------
class DNACNN(nn.Module):
    def __init__(self, in_ch=4, dropout=0.3, padding_mode="replicate"):
        super().__init__()
        self.features = nn.Sequential(
            nn.Conv1d(in_ch, 64, kernel_size=7, padding=3, padding_mode=padding_mode),
            nn.BatchNorm1d(64),
            nn.ReLU(),
            nn.MaxPool1d(2),

            nn.Conv1d(64, 128, kernel_size=7, padding=3, padding_mode=padding_mode),
            nn.BatchNorm1d(128),
            nn.ReLU(),
            nn.MaxPool1d(2),

            nn.Conv1d(128, 256, kernel_size=7, padding=6, dilation=2, padding_mode=padding_mode),
            nn.BatchNorm1d(256),
            nn.ReLU(),
        )
        self.head = nn.Sequential(
            nn.Flatten(),
            nn.Dropout(dropout),
            nn.Linear(512, 1),
        )

    def masked_pool(self, x, lens):
        B, C, Lp = x.shape
        lens_p = (lens + 3) // 4
        lens_p = torch.clamp(lens_p, min=1, max=Lp)
        t = torch.arange(Lp, device=x.device).view(1, 1, Lp)
        m = (t < lens_p.view(B, 1, 1))
        x_max = x.masked_fill(~m, -1e4).amax(dim=2, keepdim=True)
        x_sum = (x * m).sum(dim=2, keepdim=True)
        denom = lens_p.view(B, 1, 1).to(x.dtype)
        x_avg = x_sum / denom
        return x_max, x_avg

    def forward(self, x, lens):
        x = self.features(x)
        mx, av = self.masked_pool(x, lens)
        x = torch.cat([mx, av], dim=1)
        x = self.head(x)
        return x.squeeze(1)

def load_model(weights_path: str) -> nn.Module:
    model = DNACNN(dropout=0.30, padding_mode="replicate").to(device)
    sd = torch.load(weights_path, map_location=device)
    if isinstance(sd, dict) and "model_state" in sd:
        sd = sd["model_state"]
    model.load_state_dict(sd, strict=True)
    model.eval()
    return model

# -------------------------
# -------------------------
class SeqDataset(Dataset):
    def __init__(self, X_np: np.ndarray, lens_np: np.ndarray, seqs: list[str]):
        self.X = X_np
        self.lens = lens_np
        self.seqs = seqs

    def __len__(self):
        return len(self.seqs)

    def __getitem__(self, idx):
        return (
            torch.from_numpy(self.X[idx]),
            torch.as_tensor(self.lens[idx], dtype=torch.long),
            self.seqs[idx],
        )

# -------------------------
# -------------------------
_WINDOW_KERNEL_CACHE = {}

def sliding_window_sum_1d(score_1d: torch.Tensor, window: int) -> torch.Tensor:
    L = score_1d.numel()
    if L < window:
        return torch.empty(0, device=score_1d.device, dtype=score_1d.dtype)
    key = (window, score_1d.device, score_1d.dtype)
    w = _WINDOW_KERNEL_CACHE.get(key)
    if w is None:
        w = torch.ones((1, 1, window), device=score_1d.device, dtype=score_1d.dtype)
        _WINDOW_KERNEL_CACHE[key] = w
    x = score_1d.view(1, 1, L)
    y = torch.nn.functional.conv1d(x, w, stride=1, padding=0)
    return y.view(-1)

# -------------------------
# Integrated Gradients
# -------------------------
def integrated_gradients_batch(model, x, lens, steps=50, baseline_value=0.25):
    if steps <= 1:
        raise ValueError("steps must be > 1")
    model.eval()
    x = x.detach().float()
    lens = lens.detach().long()
    baseline = torch.full_like(x, float(baseline_value), dtype=torch.float32)

    alphas = torch.linspace(0.0, 1.0, steps, device=x.device, dtype=torch.float32)
    weights = torch.ones_like(alphas)
    weights[0] = 0.5
    weights[-1] = 0.5
    wsum = weights.sum()

    total_grad = torch.zeros_like(x, dtype=torch.float32)
    for a, w in zip(alphas, weights):
        x_interp = baseline + a * (x - baseline)
        x_interp.requires_grad_(True)
        logits = model(x_interp, lens)
        grads = torch.autograd.grad(outputs=logits.sum(), inputs=x_interp)[0]
        total_grad += w * grads.detach()

    avg_grad = total_grad / wsum
    return (x - baseline) * avg_grad

@torch.no_grad()
def importance_from_attr_batch(attr, lens, window=10, topk=3, edge_trim=5):
    B, _, _ = attr.shape
    score = attr.abs().sum(dim=1)  # (B,L)
    importance = np.zeros((B,), dtype=np.float32)
    best_pos = np.zeros((B,), dtype=np.int64)

    for i in range(B):
        Li = int(lens[i].item())
        s = score[i, :Li].clone()
        if edge_trim > 0 and s.numel() > 2 * edge_trim:
            s[:edge_trim] = 0
            s[-edge_trim:] = 0
        wscore = sliding_window_sum_1d(s, window=window)
        if wscore.numel() == 0:
            importance[i] = float(s.sum().item())
            best_pos[i] = 0
            continue
        k = min(int(topk), int(wscore.numel()))
        topk_vals, topk_idx = torch.topk(wscore, k=k)
        importance[i] = float(topk_vals.mean().item())
        best_pos[i] = int(topk_idx[0].item())
    return importance, best_pos

# -------------------------
# -------------------------
@torch.no_grad()
def score_probs_only(model, loader):
    probs_all = []
    for xb, lb, _ in tqdm(loader, desc="Scoring prob (AMP)", unit="batch"):
        xb = xb.to(device, non_blocking=True)
        lb = lb.to(device, non_blocking=True).view(-1)
        with autocast(amp_device, enabled=use_amp):
            p = torch.sigmoid(model(xb, lb))
        probs_all.append(p.detach().float().cpu().numpy())
    return np.concatenate(probs_all, axis=0)

# -------------------------
# -------------------------
def screen_test_by_final_score_microbatch(
    model, loader, prob_thresh, top_pct_by_final,
    window, topk_windows, ig_steps, edge_trim,
    baseline_value, ig_batch
):
    if not (0.0 < float(top_pct_by_final) <= 1.0):
        raise ValueError("top_pct_by_final must be in (0,1].")

    dev = next(model.parameters()).device
    model.eval()

    cand = []
    total = 0

    pbar1 = tqdm(loader, desc=f"Stage 1 (prob thr={prob_thresh:.2f})", unit="batch")
    for xb, lb, seqb in pbar1:
        total += len(seqb)
        xb = xb.to(dev, non_blocking=True)
        lb = lb.to(dev, non_blocking=True).view(-1)

        with torch.no_grad():
            with autocast(amp_device, enabled=use_amp):
                probs = torch.sigmoid(model(xb, lb))

        probs_cpu = probs.detach().float().cpu().numpy()
        idxs = np.where(probs_cpu > float(prob_thresh))[0]
        if idxs.size == 0:
            pbar1.set_postfix(total=total, cand=len(cand))
            continue

        xb_cpu = xb.detach().cpu()
        lb_cpu = lb.detach().cpu()
        for i in idxs.tolist():
            cand.append({
                "sequence": seqb[i],
                "prob": float(probs_cpu[i]),
                "x_cpu": xb_cpu[i],
                "l_cpu": lb_cpu[i],
            })
        pbar1.set_postfix(total=total, cand=len(cand))

    print(f"TEST total: {total} | candidates prob>{prob_thresh:.2f}: {len(cand)}")
    if not cand:
        return pd.DataFrame(columns=[
            "sequence","prob","importance","final_score",
            "best_window_start","best_window_end","best_window_seq"
        ])

    if dev.type == "cuda":
        t = torch.randn(1, device=dev, requires_grad=True)
        (t * t).sum().backward()
        torch.cuda.synchronize()

    records = []
    n = len(cand)
    pbar2 = tqdm(range(0, n, ig_batch), desc="Stage 2 (IG)", unit="chunk")

    for start in pbar2:
        chunk = cand[start:start + ig_batch]
        x_batch = torch.stack([c["x_cpu"] for c in chunk], dim=0).to(dev, non_blocking=True)
        l_batch = torch.stack([c["l_cpu"] for c in chunk], dim=0).to(dev, non_blocking=True).view(-1).long()

        # disable autocast like your dev script
        ctx = torch.cuda.amp.autocast(enabled=False) if dev.type == "cuda" else torch.amp.autocast("cpu", enabled=False)
        with ctx:
            attr = integrated_gradients_batch(
                model=model, x=x_batch, lens=l_batch,
                steps=ig_steps, baseline_value=baseline_value
            )

        imp_np, best_pos_np = importance_from_attr_batch(
            attr=attr, lens=l_batch,
            window=window, topk=topk_windows, edge_trim=edge_trim
        )

        for i, c in enumerate(chunk):
            seq = c["sequence"]
            best_pos = int(best_pos_np[i])
            best_end = best_pos + int(window)
            prob = float(c["prob"])
            importance = float(imp_np[i])
            records.append({
                "sequence": seq,
                "prob": prob,
                "importance": importance,
                "final_score": prob * importance,
                "best_window_start": best_pos,
                "best_window_end": best_end,
                "best_window_seq": seq[best_pos:best_end] if len(seq) >= best_end else "",
            })

        pbar2.set_postfix(done=min(start + ig_batch, n), total=n)

    out = pd.DataFrame(records).sort_values("final_score", ascending=False).reset_index(drop=True)
    k = max(1, int(np.ceil(len(out) * float(top_pct_by_final))))
    return out.head(k).reset_index(drop=True)

# ============================================================
# RUN
# ============================================================
if __name__ == "__main__":
    df_test = pd.read_csv(TEST_CSV)
    if SEQ_COL not in df_test.columns:
        raise ValueError(f"{TEST_CSV} missing column '{SEQ_COL}'")
    sequences = df_test[SEQ_COL].astype(str).tolist()

    print("Encoding TEST once (reuse for all reps) ...")
    X, lens = encode_sequences(sequences, seq_len=SEQ_LEN, n_fill=N_FILL)
    print("len min/med/max:", int(lens.min()), int(np.median(lens)), int(lens.max()))

    dataset = SeqDataset(X, lens, sequences)
    loader = DataLoader(dataset, batch_size=BATCH_SIZE, shuffle=False,
                        num_workers=NUM_WORKERS, pin_memory=PIN_MEMORY)

    for rep in range(1, 11):
        rep_dir = os.path.join(WEIGHTS_ROOT, f"rep{rep}")
        weights_path = os.path.join(rep_dir, "best_model.pt")
        if not os.path.exists(weights_path):
            print(f"[SKIP] rep{rep} missing weights: {weights_path}")
            continue

        thr = float(REP_THRESH[rep])
        print(f"\n==================== REP {rep} TEST (thr={thr:.2f}) ====================")

        model = load_model(weights_path)

        # --- prob only output ---
        if RUN_PROB_ONLY:
            probs = score_probs_only(model, loader)
            out_prob = pd.DataFrame({
                "sequence": sequences,
                "true_len": lens,
                "prob": probs
            })
            out_prob_csv = os.path.join(OUT_DIR, f"rep{rep}_test_probs.csv")
            out_prob.to_csv(out_prob_csv, index=False)
            print(f"Saved: {out_prob_csv}")

        # --- IG + final score output ---
        if RUN_WITH_IG:
            top_df = screen_test_by_final_score_microbatch(
                model=model,
                loader=loader,
                prob_thresh=thr,                 
                top_pct_by_final=TOP_PCT_BY_FINAL,
                window=WINDOW,
                topk_windows=TOPK_WINDOWS,
                ig_steps=IG_STEPS,
                edge_trim=EDGE_TRIM,
                baseline_value=BASELINE_VALUE,
                ig_batch=IG_BATCH,
            )
            out_ig_csv = os.path.join(OUT_DIR, f"rep{rep}_test_top_by_final.csv")
            top_df.to_csv(out_ig_csv, index=False)
            print(f"Saved: {out_ig_csv}  rows={len(top_df)}")

    print("\nAll done. Outputs in:", os.path.abspath(OUT_DIR))