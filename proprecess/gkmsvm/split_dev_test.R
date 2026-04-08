
##############
##############
library(dplyr)
all1  <- read.csv("all.csv")
all1 
stopifnot(all(c("seqnames", "label") %in% colnames(all1)))

hold_prop  <- 0.20
val_prop   <- 0.16
train_prop <- 0.64

B <- 5000
pos_label_value <- 1
n_repeat <- 10

stopifnot(abs(hold_prop + val_prop + train_prop - 1) < 1e-9)
chr_sum_full <- all1 %>%
  group_by(seqnames) %>%
  summarise(
    n = n(),
    pos = sum(label == pos_label_value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(pos_rate = pos / n)

global_pos <- mean(all1$label == pos_label_value, na.rm = TRUE)


pick_chr_set <- function(chr_sum_sub, target_n, ref_pos, B = 5000, w_pos = 2) {
  chrs <- chr_sum_sub$seqnames
  best <- list(score = Inf, sel_chr = character(0), n = 0, pos_rate = NA_real_)
  
  for (b in 1:B) {
    cand <- sample(chrs)
    
    sel <- character(0)
    cum <- 0L
    
    for (c in cand) {
      c_n <- chr_sum_sub$n[chr_sum_sub$seqnames == c]
      

      if (abs(cum - target_n) < abs((cum + c_n) - target_n)) break
      
      sel <- c(sel, c)
      cum <- cum + c_n
      if (cum >= target_n) break
    }
    
    if (length(sel) == 0) next
    
    n_sel <- sum(chr_sum_sub$n[chr_sum_sub$seqnames %in% sel])
    pos_sel <- sum(chr_sum_sub$pos[chr_sum_sub$seqnames %in% sel])
    pos_rate_sel <- pos_sel / n_sel
    
    score <- abs(n_sel - target_n) / target_n + w_pos * abs(pos_rate_sel - ref_pos)
    
    if (is.finite(score) && score < best$score) {
      best <- list(score = score, sel_chr = sel, n = n_sel, pos_rate = pos_rate_sel)
    }
  }
  
  best
}


N_total <- nrow(all1)
hold_n  <- round(N_total * hold_prop)
val_n   <- round(N_total * val_prop)


rep_chr_sets <- vector("list", n_repeat)
names(rep_chr_sets) <- paste0("rep", 1:n_repeat)

# =============================
# =============================
for (rep_i in 1:n_repeat) {
  
  set.seed(100 + rep_i)
  
  chr_sum <- chr_sum_full
  

  best_hold <- pick_chr_set(chr_sum, hold_n, global_pos, B)
  hold_chr <- best_hold$sel_chr
  

  chr_sum_remain <- chr_sum %>% filter(!(seqnames %in% hold_chr))
  remain_rows <- all1 %>% filter(!(seqnames %in% hold_chr))
  remain_pos <- mean(remain_rows$label == pos_label_value, na.rm = TRUE)
  
  best_val <- pick_chr_set(chr_sum_remain, val_n, remain_pos, B)
  val_chr <- best_val$sel_chr
  

  train_chr <- setdiff(chr_sum$seqnames, union(hold_chr, val_chr))
  

  rep_chr_sets[[rep_i]] <- list(
    hold_chr  = sort(unique(hold_chr)),
    val_chr   = sort(unique(val_chr)),
    train_chr = sort(unique(train_chr))
  )
  

  split_data <- all1 %>%
    mutate(
      set = case_when(
        seqnames %in% hold_chr ~ "hold",
        seqnames %in% val_chr  ~ "val",
        TRUE                   ~ "train"
      )
    )
  
  holdout <- split_data %>% filter(set == "hold")
  val     <- split_data %>% filter(set == "val")
  train64 <- split_data %>% filter(set == "train")
  
  # -----------------------------
  rep_dir <- paste0("rep", rep_i)
  if (!dir.exists(rep_dir)) dir.create(rep_dir, recursive = TRUE)
  

  write.csv(train64, file.path(rep_dir, "train.csv"), row.names = FALSE)
  write.csv(val,     file.path(rep_dir, "val.csv"),   row.names = FALSE)
  write.csv(holdout, file.path(rep_dir, "hold.csv"),  row.names = FALSE)
  

  cat("\n========== REP", rep_i, "==========\n")
  cat("Counts:", nrow(holdout), nrow(val), nrow(train64), "\n")
  cat("Pos rates:",
      mean(holdout$label == pos_label_value, na.rm = TRUE),
      mean(val$label == pos_label_value, na.rm = TRUE),
      mean(train64$label == pos_label_value, na.rm = TRUE), "\n")
}

# =============================
# =============================

is_same_split <- function(a, b) {
  identical(a$hold_chr,  b$hold_chr) &&
    identical(a$val_chr,   b$val_chr) &&
    identical(a$train_chr, b$train_chr)
}

dup_pairs <- list()
k <- 1
for (i in 1:(n_repeat - 1)) {
  for (j in (i + 1):n_repeat) {
    if (is_same_split(rep_chr_sets[[i]], rep_chr_sets[[j]])) {
      dup_pairs[[k]] <- c(i, j)
      k <- k + 1
    }
  }
}

cat("\n====================\n")
cat("Pairwise exact-duplicate check\n")
cat("====================\n")
if (length(dup_pairs) == 0) {
  cat("No exact duplicate splits found among reps.\n")
} else {
  cat("Exact duplicate splits found in pairs:\n")
  print(dup_pairs)
}

# =============================
# =============================
jaccard <- function(a, b) {
  if (length(union(a, b)) == 0) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}

cat("\n====================\n")
cat("Pairwise Jaccard similarity (hold / val / train)\n")
cat("====================\n")

for (i in 1:(n_repeat - 1)) {
  for (j in (i + 1):n_repeat) {
    h <- jaccard(rep_chr_sets[[i]]$hold_chr,  rep_chr_sets[[j]]$hold_chr)
    v <- jaccard(rep_chr_sets[[i]]$val_chr,   rep_chr_sets[[j]]$val_chr)
    t <- jaccard(rep_chr_sets[[i]]$train_chr, rep_chr_sets[[j]]$train_chr)
    cat(sprintf("rep%d vs rep%d: hold=%.3f  val=%.3f  train=%.3f\n", i, j, h, v, t))
  }
}

all_chr <- sort(unique(chr_sum_full$seqnames))
total_chr <- length(all_chr)

hold_union  <- sort(unique(unlist(lapply(rep_chr_sets, `[[`, "hold_chr"))))
val_union   <- sort(unique(unlist(lapply(rep_chr_sets, `[[`, "val_chr"))))
train_union <- sort(unique(unlist(lapply(rep_chr_sets, `[[`, "train_chr"))))

cat("\n====================\n")
cat("Chromosome coverage across reps\n")
cat("====================\n")
cat("Total unique chromosomes in dataset:", total_chr, "\n\n")

cat(sprintf("HOLD  coverage: %d / %d = %.1f%%\n",
            length(hold_union), total_chr, 100 * length(hold_union) / total_chr))
cat(sprintf("VAL   coverage: %d / %d = %.1f%%\n",
            length(val_union), total_chr, 100 * length(val_union) / total_chr))
cat(sprintf("TRAIN coverage: %d / %d = %.1f%%\n",
            length(train_union), total_chr, 100 * length(train_union) / total_chr))


hold_never  <- setdiff(all_chr, hold_union)
val_never   <- setdiff(all_chr, val_union)
train_never <- setdiff(all_chr, train_union)

cat("\n--- Never appeared in HOLD ---\n")
print(hold_never)
cat("\n--- Never appeared in VAL ---\n")
print(val_never)
cat("\n--- Never appeared in TRAIN ---\n")
print(train_never)

# =============================
# =============================

count_in_set <- function(chr, set_name) {
  sum(sapply(rep_chr_sets, function(x) chr %in% x[[set_name]]))
}

coverage_table <- data.frame(
  seqnames = all_chr,
  hold_count  = sapply(all_chr, count_in_set, set_name = "hold_chr"),
  val_count   = sapply(all_chr, count_in_set, set_name = "val_chr"),
  train_count = sapply(all_chr, count_in_set, set_name = "train_chr")
) %>%
  mutate(total = hold_count + val_count + train_count) %>%
  arrange(desc(hold_count), desc(val_count), desc(train_count), seqnames)

cat("\n====================\n")
cat("Per-chromosome appearance counts (0..n_repeat)\n")
cat("====================\n")
print(coverage_table)

cat("\n--- Chromosomes never in HOLD (hold_count==0) ---\n")
print(coverage_table %>% filter(hold_count == 0) %>% select(seqnames, hold_count, val_count, train_count))