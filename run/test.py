import os
import random
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from torch.amp import autocast, GradScaler
from sklearn.metrics import (
    roc_auc_score,
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    average_precision_score,
    balanced_accuracy_score,
    matthews_corrcoef,
    roc_curve,
    precision_recall_curve,
)


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
print("Device:", device)
if device.type == "cuda":
    try:
        print("GPU:", torch.cuda.get_device_name(0))
    except Exception:
        print("GPU: (unknown)")
else:
    print("GPU: None")

# -------------------------
# Config
# -------------------------
REP = 1
BASE_IN = f"../data/rep{REP}"
BASE_OUT = f"../out/rep{REP}"
os.makedirs(BASE_OUT, exist_ok=True)

TRAIN_CSV_PATH = f"{BASE_IN}/train.csv"
VAL_CSV_PATH = f"{BASE_IN}/val.csv"
HOLD_CSV_PATH = f"{BASE_IN}/hold.csv"

ROC_SAVE_PATH = f"{BASE_OUT}/roc_hold.csv"
PR_SAVE_PATH = f"{BASE_OUT}/pr_hold.csv"
HOLD_PROB_SAVE_PATH = f"{BASE_OUT}/pred_hold.csv"
best_path = f"{BASE_OUT}/best_model.pt"

SEQ_LEN = 256
BATCH_SIZE = 1024
EPOCHS = 100

# augmentation
USE_RC_AUG = False
RC_PROB = 0.5
USE_SHIFT_AUG = True
MAX_SHIFT = 5

# training controls
LR = 1e-3
WEIGHT_DECAY = 1e-4
CLIP_NORM = 1.0
DROPOUT = 0.30

# imbalance controls
POS_WEIGHT_CAP = 50.0

# early stopping
PATIENCE = 8

# DNA-friendly "N" embedding in 4ch one-hot
N_FILL = 0.25

# -------------------------
# Fast one-hot LUT
# -------------------------
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

    return x, true_len


def reverse_complement_onehot_np(x: np.ndarray) -> np.ndarray:
    xr = x[:, ::-1].copy()
    xr = xr[[3, 2, 1, 0], :]
    return xr


def random_shift_onehot_np(
    x: np.ndarray, max_shift: int, fill: float = 0.25
) -> np.ndarray:
    if max_shift <= 0:
        return x
    shift = np.random.randint(-max_shift, max_shift + 1)
    if shift == 0:
        return x

    L = x.shape[1]
    out = np.full_like(x, fill, dtype=np.float32)
    if shift > 0:
        out[:, shift:] = x[:, : L - shift]
    else:
        s = -shift
        out[:, : L - s] = x[:, s:]
    return out


# -------------------------
# Load & precompute
# -------------------------
def load_and_encode(csv_path: str, seq_len: int, n_fill: float):
    df = pd.read_csv(csv_path)
    df["sequence"] = df["sequence"].astype(str)
    df["label"] = df["label"].astype(int)

    print(f"Loaded {csv_path}: {len(df)} rows")
    print("Label value counts:\n", df["label"].value_counts(dropna=False))

    allowed = set("ACGTNacgtn")
    sample_n = min(5000, len(df))
    bad = 0
    for s in df["sequence"].head(sample_n).values:
        if any((c not in allowed) for c in s):
            bad += 1
    if sample_n > 0:
        print(
            f"Sanity check (first {sample_n} seq): {bad} contain non-ACGTN chars (treated as N-fill={n_fill})."
        )

    print(f"Precomputing one-hot + lengths for {csv_path} ...")
    tmp = [
        one_hot_encode_fast_with_len(s, seq_len, n_fill=n_fill)
        for s in df["sequence"].values
    ]
    X = np.stack([t[0] for t in tmp])
    lens = np.array([t[1] for t in tmp], dtype=np.int64)
    y = df["label"].values.astype(np.float32)
    return X, y, lens


X_train, y_train, len_train = load_and_encode(TRAIN_CSV_PATH, SEQ_LEN, N_FILL)
X_val, y_val, len_val = load_and_encode(VAL_CSV_PATH, SEQ_LEN, N_FILL)

print("Train/Val:", len(y_train), len(y_val))
print("pos rate train:", float(y_train.mean()), "val:", float(y_val.mean()))
print(
    "len train min/med/max:",
    int(len_train.min()),
    int(np.median(len_train)),
    int(len_train.max()),
)
print(
    "len val   min/med/max:",
    int(len_val.min()),
    int(np.median(len_val)),
    int(len_val.max()),
)


class DNADataset(Dataset):
    def __init__(
        self,
        X_np,
        y_np,
        len_np,
        train,
        use_rc,
        rc_prob,
        use_shift,
        max_shift,
        n_fill=0.25,
    ):
        self.X = X_np
        self.y = y_np
        self.len = len_np
        self.train = train
        self.use_rc = use_rc
        self.rc_prob = rc_prob
        self.use_shift = use_shift
        self.max_shift = max_shift
        self.n_fill = n_fill

    def __len__(self):
        return self.y.shape[0]

    def __getitem__(self, idx):
        x = self.X[idx]
        y = self.y[idx]
        Ltrue = int(self.len[idx])

        if self.train:
            if self.use_rc and (np.random.rand() < self.rc_prob):
                x = reverse_complement_onehot_np(x)
            if self.use_shift and self.max_shift > 0:
                x = random_shift_onehot_np(x, self.max_shift, fill=self.n_fill)

        return (
            torch.from_numpy(x),
            torch.as_tensor(y, dtype=torch.float32),
            torch.as_tensor(Ltrue, dtype=torch.long),
        )


# -------------------------
# DataLoaders
# -------------------------
num_workers = 8 if os.name != "nt" else 0
pin = device.type == "cuda"

dl_kwargs = dict(batch_size=BATCH_SIZE, pin_memory=pin)
if num_workers > 0:
    dl_kwargs.update(
        dict(num_workers=num_workers, persistent_workers=True, prefetch_factor=4)
    )
else:
    dl_kwargs.update(dict(num_workers=0))

train_loader = DataLoader(
    DNADataset(
        X_train,
        y_train,
        len_train,
        train=True,
        use_rc=USE_RC_AUG,
        rc_prob=RC_PROB,
        use_shift=USE_SHIFT_AUG,
        max_shift=MAX_SHIFT,
        n_fill=N_FILL,
    ),
    shuffle=True,
    **dl_kwargs,
)

val_loader = DataLoader(
    DNADataset(
        X_val,
        y_val,
        len_val,
        train=False,
        use_rc=False,
        rc_prob=0.0,
        use_shift=False,
        max_shift=0,
        n_fill=N_FILL,
    ),
    shuffle=False,
    **dl_kwargs,
)


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
            nn.Conv1d(
                128,
                256,
                kernel_size=7,
                padding=6,
                dilation=2,
                padding_mode=padding_mode,
            ),
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
        m = t < lens_p.view(B, 1, 1)

        x_max = x.masked_fill(~m, float("-inf")).amax(dim=2, keepdim=True)
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


model = DNACNN(dropout=DROPOUT, padding_mode="replicate").to(device)
optimizer = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

pos = float(y_train.sum())
neg = float(len(y_train) - y_train.sum())
pw = neg / max(pos, 1.0)
pw = min(pw, POS_WEIGHT_CAP)
pos_weight = torch.tensor([pw], device=device, dtype=torch.float32)
print(f"pos_weight used: {float(pos_weight.item()):.4f}")

criterion = nn.BCEWithLogitsLoss(pos_weight=pos_weight)

use_amp = device.type == "cuda"
amp_device = "cuda" if use_amp else "cpu"
scaler = GradScaler(amp_device, enabled=use_amp)

scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
    optimizer, mode="max", factor=0.5, patience=2
)


def safe_auc(y_true: np.ndarray, y_prob: np.ndarray) -> float:
    if len(np.unique(y_true)) < 2:
        return float("nan")
    return roc_auc_score(y_true, y_prob)


def best_threshold_by_metric(y_true: np.ndarray, y_prob: np.ndarray, metric="f1"):
    thresholds = np.linspace(0.0, 1.0, 101)
    best_t, best_v = 0.5, -1.0
    for t in thresholds:
        pred = (y_prob >= t).astype(np.int32)
        if metric == "f1":
            v = f1_score(y_true, pred, zero_division=0)
        elif metric == "bal_acc":
            v = balanced_accuracy_score(y_true, pred)
        else:
            v = accuracy_score(y_true, pred)
        if v > best_v:
            best_v, best_t = v, t
    return best_t, best_v


@torch.no_grad()
def evaluate(model: nn.Module, loader: DataLoader):
    model.eval()
    ys, ps = [], []
    for x, yb, lb in loader:
        x = x.to(device, non_blocking=True)
        yb = yb.to(device, non_blocking=True)
        lb = lb.to(device, non_blocking=True)

        logits = model(x, lb)
        probs = torch.sigmoid(logits)

        ys.append(yb.detach().cpu())
        ps.append(probs.detach().cpu())

    y_true = torch.cat(ys).numpy().astype(np.int32)
    y_prob = torch.cat(ps).numpy()

    auc = safe_auc(y_true, y_prob)
    ap = average_precision_score(y_true, y_prob)

    def cls_metrics(threshold: float):
        pred = (y_prob >= threshold).astype(np.int32)
        return {
            "acc": float(accuracy_score(y_true, pred)),
            "precision": float(precision_score(y_true, pred, zero_division=0)),
            "recall": float(recall_score(y_true, pred, zero_division=0)),
            "f1": float(f1_score(y_true, pred, zero_division=0)),
            "bal_acc": float(balanced_accuracy_score(y_true, pred)),
            "mcc": float(matthews_corrcoef(y_true, pred)),
        }

    m05 = cls_metrics(0.5)
    t_f1, best_f1 = best_threshold_by_metric(y_true, y_prob, metric="f1")
    mf1 = cls_metrics(t_f1)
    t_bal, best_bal = best_threshold_by_metric(y_true, y_prob, metric="bal_acc")
    mbal = cls_metrics(t_bal)

    q = np.quantile(y_prob, [0, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0]).astype(float)

    return {
        "auc": float(auc),
        "ap": float(ap),
        "acc@0.5": m05["acc"],
        "prec@0.5": m05["precision"],
        "recall@0.5": m05["recall"],
        "f1@0.5": m05["f1"],
        "bal@0.5": m05["bal_acc"],
        "mcc@0.5": m05["mcc"],
        "t_f1": float(t_f1),
        "best_f1": float(best_f1),
        "acc@t_f1": mf1["acc"],
        "prec@t_f1": mf1["precision"],
        "recall@t_f1": mf1["recall"],
        "f1@t_f1": mf1["f1"],
        "bal@t_f1": mf1["bal_acc"],
        "mcc@t_f1": mf1["mcc"],
        "t_bal": float(t_bal),
        "best_bal": float(best_bal),
        "acc@t_bal": mbal["acc"],
        "prec@t_bal": mbal["precision"],
        "recall@t_bal": mbal["recall"],
        "f1@t_bal": mbal["f1"],
        "p_q": q,
    }


# -------------------------
# Train
# -------------------------
best_auc = -1.0
no_improve = 0

for epoch in range(1, EPOCHS + 1):
    model.train()
    running_loss = 0.0
    n_batches = 0

    for x, yb, lb in train_loader:
        x = x.to(device, non_blocking=True)
        yb = yb.to(device, non_blocking=True)
        lb = lb.to(device, non_blocking=True)

        optimizer.zero_grad(set_to_none=True)

        with autocast(amp_device, enabled=use_amp):
            logits = model(x, lb)
            loss = criterion(logits, yb)

        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.parameters(), CLIP_NORM)
        scaler.step(optimizer)
        scaler.update()

        running_loss += float(loss.detach().item())
        n_batches += 1

    avg_loss = running_loss / max(n_batches, 1)
    metrics = evaluate(model, val_loader)
    auc = metrics["auc"]

    if not np.isnan(auc):
        scheduler.step(auc)

    q = metrics["p_q"]

    improved = (not np.isnan(auc)) and (auc > best_auc)
    if improved:
        best_auc = auc
        no_improve = 0
        torch.save(model.state_dict(), best_path)
        print(
            f"BEST @ Epoch {epoch:02d} | loss={avg_loss:.4f} | "
            f"AUC={metrics['auc']:.4f} | AP={metrics['ap']:.4f} | "
            f"@0.5 acc={metrics['acc@0.5']:.4f} prec={metrics['prec@0.5']:.4f} recall={metrics['recall@0.5']:.4f} f1={metrics['f1@0.5']:.4f} | "
            f"@t_f1({metrics['t_f1']:.2f}) acc={metrics['acc@t_f1']:.4f} prec={metrics['prec@t_f1']:.4f} recall={metrics['recall@t_f1']:.4f} f1={metrics['f1@t_f1']:.4f} | "
            f"bal@0.5={metrics['bal@0.5']:.4f} mcc@0.5={metrics['mcc@0.5']:.4f} | "
            f"best_bal={metrics['best_bal']:.4f} (t_bal={metrics['t_bal']:.2f}) | "
            f"p_q=[{q[0]:.3f},{q[1]:.3f},{q[2]:.3f},{q[3]:.3f},{q[4]:.3f},{q[5]:.3f},{q[6]:.3f}]"
        )
    else:
        if not np.isnan(auc):
            no_improve += 1

    if no_improve >= PATIENCE:
        print(
            f"Early stopping: no AUC improvement for {PATIENCE} epochs. Best AUC={best_auc:.4f}"
        )
        break


# -------------------------
# Holdout evaluation (threshold picked by F1 on VAL)
# -------------------------
X_hold, y_hold, len_hold = load_and_encode(HOLD_CSV_PATH, SEQ_LEN, N_FILL)

hold_loader = DataLoader(
    DNADataset(
        X_hold,
        y_hold,
        len_hold,
        train=False,
        use_rc=False,
        rc_prob=0.0,
        use_shift=False,
        max_shift=0,
        n_fill=N_FILL,
    ),
    shuffle=False,
    **dl_kwargs,
)

assert os.path.exists(best_path), f"Missing {best_path}. Train first or check path."
model.load_state_dict(torch.load(best_path, map_location=device))
model.eval()
print(f"\nLoaded best weights from: {best_path}")


@torch.no_grad()
def predict_probs(model: nn.Module, loader: DataLoader):
    model.eval()
    ys, ps, ls = [], [], []
    for x, yb, lb in loader:
        x = x.to(device, non_blocking=True)
        lb = lb.to(device, non_blocking=True)
        logits = model(x, lb)
        probs = torch.sigmoid(logits)
        ys.append(yb.detach().cpu())
        ps.append(probs.detach().cpu())
        ls.append(lb.detach().cpu())
    y_true = torch.cat(ys).numpy().astype(np.int32)
    y_prob = torch.cat(ps).numpy().astype(np.float64)
    lens = torch.cat(ls).numpy().astype(np.int64)
    return y_true, y_prob, lens


y_val_true, y_val_prob, _ = predict_probs(model, val_loader)
t_f1, best_f1 = best_threshold_by_metric(y_val_true, y_val_prob, metric="f1")
print(f"\n[F1-threshold on VAL] t_f1={t_f1:.4f}, best_f1(val)={best_f1:.4f}")

y_hold_true, y_hold_prob, hold_lens = predict_probs(model, hold_loader)

hold_auc = safe_auc(y_hold_true, y_hold_prob)
hold_ap = average_precision_score(y_hold_true, y_hold_prob)

hold_pred = (y_hold_prob >= t_f1).astype(np.int32)

hold_metrics = {
    "auc": float(hold_auc),
    "ap": float(hold_ap),
    "threshold": float(t_f1),
    "acc": float(accuracy_score(y_hold_true, hold_pred)),
    "precision": float(precision_score(y_hold_true, hold_pred, zero_division=0)),
    "recall": float(recall_score(y_hold_true, hold_pred, zero_division=0)),
    "f1": float(f1_score(y_hold_true, hold_pred, zero_division=0)),
    "bal_acc": float(balanced_accuracy_score(y_hold_true, hold_pred)),
    "mcc": float(matthews_corrcoef(y_hold_true, hold_pred)),
}

print("\n[HOLD metrics using VAL F1 threshold]")
for k, v in hold_metrics.items():
    print(f"{k:>10s}: {v:.6f}" if isinstance(v, float) else f"{k:>10s}: {v}")

fpr, tpr, roc_th = roc_curve(y_hold_true, y_hold_prob)
roc_df = pd.DataFrame({"fpr": fpr, "tpr": tpr, "threshold": roc_th})
roc_df.to_csv(ROC_SAVE_PATH, index=False)
print(f"\nSaved ROC data to: {ROC_SAVE_PATH}  (columns: fpr,tpr,threshold)")

prec, rec, pr_th = precision_recall_curve(y_hold_true, y_hold_prob)
pr_df = pd.DataFrame(
    {
        "precision": prec,
        "recall": rec,
        "threshold": np.r_[pr_th, np.nan],
    }
)
pr_df.to_csv(PR_SAVE_PATH, index=False)
print(f"Saved PR data to:  {PR_SAVE_PATH}  (columns: precision,recall,threshold)")

pred_df = pd.DataFrame(
    {
        "y_true": y_hold_true.astype(int),
        "y_prob": y_hold_prob.astype(float),
        "y_pred": hold_pred.astype(int),
        "len_true": hold_lens.astype(int),
    }
)
pred_df.to_csv(HOLD_PROB_SAVE_PATH, index=False)
print(f"Saved per-sample preds to: {HOLD_PROB_SAVE_PATH}")
