# ==========================================================
# gkmSVM end-to-end in R (BY-REP, model selection on VAL AUC)
# - Loop rep1..repN (train/val/hold)
# - Hyperparameter grid search within each rep
# - Model selected ONLY on VAL AUC
# - Threshold picked ONLY on VAL by best F1
# - Evaluate ONLY on HOLD
# - Save per-rep artifacts:
#     hold_roc_curve.png
#     hold_roc_data.csv
#     hold_confusion_matrix.csv
#     hold_report.txt
#     hold_summary.json
#     val_scores_best.txt / hold_scores_best.txt
#     candidate_models_summary.csv
# - Save global summaries:
#     all_reps_summary.csv
#     all_reps_roc_summary.csv
#     mean_roc_curve.csv
#     all_reps_auc_summary.csv
# ==========================================================

suppressPackageStartupMessages({
  library(gkmSVM)
  library(data.table)
  library(pROC)
  library(PRROC)
})

SEED <- 123
set.seed(SEED)

# -----------------------------
# CONFIG
# -----------------------------
REPS_ROOT <- "data"                 # contains rep1..repN
OUT_ROOT  <- "./gkmsvm_by_rep_grid" # output root
N_REP     <- 10

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# HYPERPARAMETER GRID
# Adjust based on your gkmSVM version / runtime budget
# -----------------------------
PARAM_GRID <- list(
  list(L = 10, K = 6, maxnmm = 3),
  list(L = 10, K = 5, maxnmm = 3),
  list(L = 10, K = 6, maxnmm = 2),
  list(L = 8,  K = 5, maxnmm = 2),
  list(L = 8,  K = 6, maxnmm = 2)
)

# -----------------------------
# HELPERS
# -----------------------------
clean_seq <- function(s) {
  s <- toupper(as.character(s))
  s <- gsub("[^ACGT]", "N", s)
  s
}

write_fasta <- function(seqs, file, prefix) {
  con <- file(file, open = "w")
  on.exit(close(con), add = TRUE)
  for (i in seq_along(seqs)) {
    writeLines(paste0(">", prefix, "_", i), con)
    writeLines(seqs[[i]], con)
  }
}

read_score_file <- function(score_file) {
  dt <- data.table::fread(score_file, header = FALSE)
  as.numeric(dt[[ncol(dt)]])
}

# pick threshold on VAL by max F1
best_thresh_f1 <- function(y, score, grid) {
  best_t <- median(score, na.rm = TRUE)
  best_f1 <- -Inf
  best_prec <- NA_real_
  best_rec  <- NA_real_
  
  for (t in grid) {
    pred <- ifelse(score >= t, 1L, 0L)
    tp <- sum(y == 1L & pred == 1L)
    fp <- sum(y == 0L & pred == 1L)
    fn <- sum(y == 1L & pred == 0L)
    
    prec <- tp / (tp + fp + 1e-12)
    rec  <- tp / (tp + fn + 1e-12)
    f1 <- 2 * prec * rec / (prec + rec + 1e-12)
    
    if (f1 > best_f1) {
      best_f1 <- f1
      best_t <- t
      best_prec <- prec
      best_rec <- rec
    }
  }
  
  list(
    threshold = best_t,
    f1 = best_f1,
    precision = best_prec,
    recall = best_rec
  )
}

# mean ROC (interpolation) for multiple reps
mean_roc_from_list <- function(roc_list, fpr_grid = seq(0, 1, length.out = 200)) {
  tpr_mat <- matrix(NA_real_, nrow = length(roc_list), ncol = length(fpr_grid))
  
  for (i in seq_along(roc_list)) {
    df <- roc_list[[i]]
    df <- df[order(df$fpr), ]
    
    fpr <- pmin(pmax(df$fpr, 0), 1)
    tpr <- pmin(pmax(df$tpr, 0), 1)
    
    # remove duplicated fpr values for approx
    keep <- !duplicated(fpr)
    fpr <- fpr[keep]
    tpr <- tpr[keep]
    
    if (length(fpr) < 2) {
      next
    }
    
    tpr_i <- approx(x = fpr, y = tpr, xout = fpr_grid, rule = 2)$y
    tpr_i[1] <- 0
    tpr_mat[i, ] <- tpr_i
  }
  
  mean_tpr <- colMeans(tpr_mat, na.rm = TRUE)
  std_tpr  <- apply(tpr_mat, 2, sd, na.rm = TRUE)
  mean_auc <- sum(diff(fpr_grid) * (head(mean_tpr, -1) + tail(mean_tpr, -1)) / 2)
  
  list(
    df = data.frame(fpr = fpr_grid, mean_tpr = mean_tpr, std_tpr = std_tpr),
    mean_auc = mean_auc
  )
}

safe_auc <- function(y, score) {
  roc_obj <- pROC::roc(response = y, predictor = score, quiet = TRUE, direction = "<")
  as.numeric(pROC::auc(roc_obj))
}

# -----------------------------
# MAIN: run one rep
# -----------------------------
run_one_rep <- function(rep_i) {
  rep_dir <- file.path(REPS_ROOT, paste0("rep", rep_i))
  train_csv <- file.path(rep_dir, "train.csv")
  val_csv   <- file.path(rep_dir, "val.csv")
  hold_csv  <- file.path(rep_dir, "hold.csv")
  
  if (!file.exists(train_csv) || !file.exists(val_csv) || !file.exists(hold_csv)) {
    cat(sprintf("[SKIP] rep%d missing files.\n", rep_i))
    return(NULL)
  }
  
  cat(sprintf("\n==================== REP %d ====================\n", rep_i))
  
  rep_out <- file.path(OUT_ROOT, paste0("rep", rep_i))
  dir.create(rep_out, showWarnings = FALSE, recursive = TRUE)
  
  # ---- read ----
  train <- fread(train_csv)
  val   <- fread(val_csv)
  hold  <- fread(hold_csv)
  
  stopifnot(all(c("sequence", "label") %in% names(train)))
  stopifnot(all(c("sequence", "label") %in% names(val)))
  stopifnot(all(c("sequence", "label") %in% names(hold)))
  
  train[, sequence := clean_seq(sequence)]
  val[,   sequence := clean_seq(sequence)]
  hold[,  sequence := clean_seq(sequence)]
  
  train[, label := as.integer(label)]
  val[,   label := as.integer(label)]
  hold[,  label := as.integer(label)]
  
  stopifnot(all(train$label %in% c(0L, 1L)))
  stopifnot(all(val$label   %in% c(0L, 1L)))
  stopifnot(all(hold$label  %in% c(0L, 1L)))
  
  y_val  <- val$label
  y_hold <- hold$label
  
  cat("Counts train/val/hold:", nrow(train), nrow(val), nrow(hold), "\n")
  cat("Pos rate train/val/hold:",
      mean(train$label == 1L), mean(val$label == 1L), mean(hold$label == 1L), "\n")
  
  # ---- write FASTA once per rep ----
  train_pos_fa <- file.path(rep_out, "train_pos.fa")
  train_neg_fa <- file.path(rep_out, "train_neg.fa")
  val_fa       <- file.path(rep_out, "val.fa")
  hold_fa      <- file.path(rep_out, "hold.fa")
  
  write_fasta(train[label == 1L, sequence], train_pos_fa, "pos")
  write_fasta(train[label == 0L, sequence], train_neg_fa, "neg")
  write_fasta(val$sequence,  val_fa,  "val")
  write_fasta(hold$sequence, hold_fa, "hold")
  
  # ---- grid search on VAL AUC ----
  candidate_rows <- list()
  best_model <- NULL
  best_val_auc <- -Inf
  
  for (j in seq_along(PARAM_GRID)) {
    p <- PARAM_GRID[[j]]
    tag <- sprintf("L%s_K%s_d%s", p$L, p$K, p$maxnmm)
    cat(sprintf("  [rep%d] Candidate %d/%d: %s\n", rep_i, j, length(PARAM_GRID), tag))
    
    kernel_fn    <- file.path(rep_out, paste0("train_kernel_", tag, ".txt"))
    model_prfx   <- file.path(rep_out, paste0("gkmsvm_model_", tag))
    val_scores_fn <- file.path(rep_out, paste0("val_scores_", tag, ".txt"))
    
    ok <- TRUE
    err_msg <- NA_character_
    
    tryCatch({
      # NOTE:
      # If your installed gkmSVM version uses different arg names,
      # change L/K/maxnmm here accordingly.
      gkmsvm_kernel(
        train_pos_fa,
        train_neg_fa,
        kernel_fn,
        L = p$L,
        K = p$K,
        maxnmm = p$maxnmm
      )
      
      gkmsvm_train(
        kernel_fn,
        train_pos_fa,
        train_neg_fa,
        model_prfx
      )
      
      gkmsvm_classify(
        val_fa,
        model_prfx,
        val_scores_fn
      )
    }, error = function(e) {
      ok <<- FALSE
      err_msg <<- conditionMessage(e)
    })
    
    if (!ok || !file.exists(val_scores_fn)) {
      candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
        rep = rep_i,
        candidate_tag = tag,
        L = p$L,
        K = p$K,
        maxnmm = p$maxnmm,
        val_auc = NA_real_,
        status = "failed",
        error_message = ifelse(is.na(err_msg), "", err_msg),
        stringsAsFactors = FALSE
      )
      next
    }
    
    val_scores <- read_score_file(val_scores_fn)
    if (length(val_scores) != nrow(val)) {
      candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
        rep = rep_i,
        candidate_tag = tag,
        L = p$L,
        K = p$K,
        maxnmm = p$maxnmm,
        val_auc = NA_real_,
        status = "failed_length_mismatch",
        error_message = sprintf("Expected %d val scores, got %d", nrow(val), length(val_scores)),
        stringsAsFactors = FALSE
      )
      next
    }
    
    val_auc <- safe_auc(y_val, val_scores)
    
    candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
      rep = rep_i,
      candidate_tag = tag,
      L = p$L,
      K = p$K,
      maxnmm = p$maxnmm,
      val_auc = val_auc,
      status = "ok",
      error_message = "",
      stringsAsFactors = FALSE
    )
    
    if (is.finite(val_auc) && val_auc > best_val_auc) {
      best_val_auc <- val_auc
      best_model <- list(
        tag = tag,
        params = p,
        model_prfx = model_prfx,
        val_scores_fn = val_scores_fn,
        val_scores = val_scores
      )
    }
  }
  
  candidate_df <- rbindlist(candidate_rows, fill = TRUE)
  candidate_df <- candidate_df[order(-val_auc)]
  candidate_csv <- file.path(rep_out, "candidate_models_summary.csv")
  fwrite(candidate_df, candidate_csv)
  
  if (is.null(best_model)) {
    cat(sprintf("[FAIL] rep%d no valid candidate model.\n", rep_i))
    return(NULL)
  }
  
  cat(sprintf("  [rep%d] Best model by VAL AUC: %s (AUC=%.4f)\n",
              rep_i, best_model$tag, best_val_auc))
  
  # ---- threshold on VAL by best F1 (keep same logic as your original request) ----
  val_scores_best <- best_model$val_scores
  val_scores_best_fn <- file.path(rep_out, "val_scores_best.txt")
  file.copy(best_model$val_scores_fn, val_scores_best_fn, overwrite = TRUE)
  
  grid_q <- as.numeric(quantile(
    val_scores_best,
    probs = seq(0.05, 0.95, by = 0.05),
    na.rm = TRUE
  ))
  grid_q <- sort(unique(grid_q))
  
  bt <- best_thresh_f1(y_val, val_scores_best, grid = grid_q)
  thresh <- bt$threshold
  
  # ---- score HOLD with best model ----
  hold_scores_best_fn <- file.path(rep_out, "hold_scores_best.txt")
  gkmsvm_classify(hold_fa, best_model$model_prfx, hold_scores_best_fn)
  
  hold_scores <- read_score_file(hold_scores_best_fn)
  stopifnot(length(hold_scores) == nrow(hold))
  
  # ---- HOLD ROC + ROC data ----
  roc_obj <- pROC::roc(response = y_hold, predictor = hold_scores, quiet = TRUE, direction = "<")
  hold_auc <- as.numeric(pROC::auc(roc_obj))
  
  roc_df <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities,
    threshold = roc_obj$thresholds
  )
  
  roc_csv <- file.path(rep_out, "hold_roc_data.csv")
  write.csv(roc_df, roc_csv, row.names = FALSE)
  
  roc_png <- file.path(rep_out, "hold_roc_curve.png")
  png(roc_png, width = 900, height = 700)
  plot(roc_obj, main = sprintf("HOLD ROC (best gkmSVM) AUC=%.4f", hold_auc))
  abline(a = 0, b = 1, lty = 2)
  dev.off()
  
  # ---- HOLD AP ----
  scores_pos <- hold_scores[y_hold == 1L]
  scores_neg <- hold_scores[y_hold == 0L]
  pr <- PRROC::pr.curve(scores.class0 = scores_pos, scores.class1 = scores_neg, curve = FALSE)
  hold_ap <- as.numeric(pr$auc.integral)
  
  # ---- HOLD confusion matrix / ACC / report (threshold from VAL F1) ----
  pred_hold <- ifelse(hold_scores >= thresh, 1L, 0L)
  hold_acc <- mean(pred_hold == y_hold)
  
  cm <- table(
    true = factor(y_hold, levels = c(0, 1)),
    pred = factor(pred_hold, levels = c(0, 1))
  )
  
  cm_csv <- file.path(rep_out, "hold_confusion_matrix.csv")
  write.csv(as.data.frame.matrix(cm), cm_csv)
  
  tp <- cm["1", "1"]
  tn <- cm["0", "0"]
  fp <- cm["0", "1"]
  fn <- cm["1", "0"]
  
  prec1 <- tp / (tp + fp + 1e-12)
  rec1  <- tp / (tp + fn + 1e-12)
  f11   <- 2 * prec1 * rec1 / (prec1 + rec1 + 1e-12)
  
  prec0 <- tn / (tn + fn + 1e-12)
  rec0  <- tn / (tn + fp + 1e-12)
  f10   <- 2 * prec0 * rec0 / (prec0 + rec0 + 1e-12)
  
  report_txt <- file.path(rep_out, "hold_report.txt")
  cat(
    sprintf("gkmSVM HOLD RESULTS (rep%d)\n", rep_i),
    "====================================\n\n",
    sprintf("Best model tag (selected on VAL AUC): %s\n", best_model$tag),
    sprintf("Best model params: L=%s, K=%s, maxnmm=%s\n",
            best_model$params$L, best_model$params$K, best_model$params$maxnmm),
    sprintf("VAL AUC (best model): %.6f\n", best_val_auc),
    sprintf("Threshold (picked on VAL by best F1): %.6f\n", thresh),
    sprintf("VAL best F1: %.6f\n", bt$f1),
    sprintf("VAL precision at best F1 threshold: %.6f\n", bt$precision),
    sprintf("VAL recall at best F1 threshold: %.6f\n\n", bt$recall),
    sprintf("HOLD AUC: %.6f\n", hold_auc),
    sprintf("HOLD AP : %.6f\n", hold_ap),
    sprintf("HOLD ACC: %.6f\n\n", hold_acc),
    "Confusion Matrix (rows=true, cols=pred)\n",
    paste(capture.output(print(cm)), collapse = "\n"),
    "\n\nPer-class metrics\n",
    sprintf("Class 1: precision=%.4f recall=%.4f f1=%.4f\n", prec1, rec1, f11),
    sprintf("Class 0: precision=%.4f recall=%.4f f1=%.4f\n", prec0, rec0, f10),
    file = report_txt,
    sep = ""
  )
  
  # ---- summary json ----
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite")
  }
  
  summary <- list(
    rep = rep_i,
    best_model_tag = best_model$tag,
    best_model_params = best_model$params,
    val_auc_best_model = best_val_auc,
    threshold_from_val_f1 = thresh,
    val_best_f1 = bt$f1,
    hold_metrics = list(
      auc = hold_auc,
      ap = hold_ap,
      acc = hold_acc,
      class1_precision = prec1,
      class1_recall = rec1,
      class1_f1 = f11,
      class0_precision = prec0,
      class0_recall = rec0,
      class0_f1 = f10
    ),
    artifacts = list(
      candidate_models_summary_csv = candidate_csv,
      hold_roc_curve_png = roc_png,
      hold_roc_data_csv = roc_csv,
      hold_confusion_matrix_csv = cm_csv,
      val_scores_best_txt = val_scores_best_fn,
      hold_scores_best_txt = hold_scores_best_fn,
      hold_report_txt = report_txt
    )
  )
  
  json_path <- file.path(rep_out, "hold_summary.json")
  jsonlite::write_json(summary, path = json_path, pretty = TRUE, auto_unbox = TRUE)
  
  cat(sprintf("REP %d saved to: %s\n", rep_i, normalizePath(rep_out)))
  cat(sprintf(
    "Best=%s | VAL AUC=%.4f | HOLD AUC=%.4f | AP=%.4f | ACC=%.4f | thr=%.6f\n",
    best_model$tag, best_val_auc, hold_auc, hold_ap, hold_acc, thresh
  ))
  
  # row for global summary
  list(
    rep = rep_i,
    best_model_tag = best_model$tag,
    best_model_L = best_model$params$L,
    best_model_K = best_model$params$K,
    best_model_maxnmm = best_model$params$maxnmm,
    val_auc_best_model = best_val_auc,
    threshold_from_val_f1 = thresh,
    val_best_f1 = bt$f1,
    hold_auc = hold_auc,
    hold_ap = hold_ap,
    hold_acc = hold_acc,
    n_hold = nrow(hold),
    pos_rate_hold = mean(y_hold == 1L),
    roc_df = roc_df
  )
}

# -----------------------------
# RUN ALL REPS + GLOBAL SUMMARIES
# -----------------------------
rows <- list()
roc_list <- list()
auc_rows <- list()

for (rep_i in 1:N_REP) {
  out <- run_one_rep(rep_i)
  if (!is.null(out)) {
    rows[[length(rows) + 1]] <- data.frame(
      rep = out$rep,
      best_model_tag = out$best_model_tag,
      best_model_L = out$best_model_L,
      best_model_K = out$best_model_K,
      best_model_maxnmm = out$best_model_maxnmm,
      val_auc_best_model = out$val_auc_best_model,
      threshold_from_val_f1 = out$threshold_from_val_f1,
      val_best_f1 = out$val_best_f1,
      hold_auc = out$hold_auc,
      hold_ap = out$hold_ap,
      hold_acc = out$hold_acc,
      n_hold = out$n_hold,
      pos_rate_hold = out$pos_rate_hold
    )
    
    roc_list[[paste0("rep", out$rep)]] <- out$roc_df
    
    auc_rows[[length(auc_rows) + 1]] <- data.frame(
      rep = out$rep,
      hold_auc = out$hold_auc
    )
  }
}

if (length(rows) == 0) {
  cat("No reps processed.\n")
  quit(save = "no")
}

# ---- all reps summary ----
summary_df <- rbindlist(rows, fill = TRUE)
summary_df <- summary_df[order(summary_df$rep), ]
summary_csv <- file.path(OUT_ROOT, "all_reps_summary.csv")
fwrite(summary_df, summary_csv)
cat("\nSaved:", normalizePath(summary_csv), "\n")

# ---- all reps roc long table ----
roc_long <- rbindlist(lapply(names(roc_list), function(nm) {
  rep_id <- as.integer(gsub("^rep", "", nm))
  df <- roc_list[[nm]]
  df$rep <- rep_id
  df
}), fill = TRUE)

roc_long <- roc_long[order(roc_long$rep, roc_long$fpr), ]
roc_long_csv <- file.path(OUT_ROOT, "all_reps_roc_summary.csv")
fwrite(roc_long, roc_long_csv)
cat("Saved:", normalizePath(roc_long_csv), "\n")

# ---- mean ROC ----
mr <- mean_roc_from_list(roc_list, fpr_grid = seq(0, 1, length.out = 200))
mean_roc_csv <- file.path(OUT_ROOT, "mean_roc_curve.csv")
fwrite(mr$df, mean_roc_csv)
cat("Saved:", normalizePath(mean_roc_csv), "\n")
cat(sprintf("Mean AUC (interpolated): %.4f\n", mr$mean_auc))

# ---- auc summary ----
auc_df <- rbindlist(auc_rows, fill = TRUE)
auc_df <- auc_df[order(auc_df$rep), ]
auc_csv <- file.path(OUT_ROOT, "all_reps_auc_summary.csv")
fwrite(auc_df, auc_csv)
cat("Saved:", normalizePath(auc_csv), "\n")

cat("\nAll artifacts saved under:", normalizePath(OUT_ROOT), "\n")
