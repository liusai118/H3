library(data.table)
process_sequence_score <- function(score_file, label_file, output_file = NULL) {
  score_dt <- fread(score_file)
  label_dt <- fread(label_file)
  required_score_cols <- c("sequence", "importance")
  required_label_cols <- c("sequence", "label_sum")
  if (!all(required_score_cols %in% names(score_dt))) {
    stop(sprintf(
      "score_file : %s",
      paste(setdiff(required_score_cols, names(score_dt)), collapse = ", ")
    ))
  }
  if (!all(required_label_cols %in% names(label_dt))) {
    stop(sprintf(
      "label_file : %s",
      paste(setdiff(required_label_cols, names(label_dt)), collapse = ", ")
    ))
  }
  score_dt <- unique(score_dt, by = "sequence")
  label_dt <- unique(label_dt, by = "sequence")
  merged_dt <- merge(
    score_dt,
    label_dt[, .(sequence, label_sum)],
    by = "sequence",
    all.x = TRUE
  )
  merged_dt[, score := label_sum * importance]
  setorder(merged_dt, -score)
  result_dt <- merged_dt[, .(sequence, importance, label_sum, score)]
  if (!is.null(output_file)) {
    fwrite(result_dt, output_file)
  }
  return(result_dt)
}

neg_result <- process_sequence_score(
  score_file = "./neg_score/average_prob_importance_per_sequence_neg.csv",
  label_file = "./test_ig_by_rep_neg/average_prob_importance_per_sequence.csv",
  output_file = "../Dataset/neg.csv"
)

pos_result <- process_sequence_score(
  score_file = "./pos_score/average_prob_importance_per_sequence_pos.csv",
  label_file = "./test_ig_by_rep_pos/average_prob_importance_per_sequence.csv",
  output_file = "../Dataset/pos.csv"
)

mut_result <- process_sequence_score(
  score_file = "./mut_score/average_prob_importance_per_sequence_mut.csv",
  label_file = "./test_ig_by_rep_mut/average_prob_importance_per_sequence.csv",
  output_file = "../Dataset/mut.csv"
)

