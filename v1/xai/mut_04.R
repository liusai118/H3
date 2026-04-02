pos <- read.csv(Dataset/pos.csv)
pos <- pos[order(pos$score,decreasing = T),]
top_sequence <- pos[1:10000,]
mean(top_sequence$score)
ggplot(top_sequence, aes(x = score)) +
  geom_histogram(bins =30, fill = "skyblue", color = "black") +
  labs(title = "Histogram of score",
       x = "score",
       y = "Frequency")
###########################################################################################
set.seed(123)
mutate_seq <- function(seq, n_mut) {
  bases <- c("A","T","C","G")
  seq_vec <- strsplit(seq, "")[[1]]
  
  pos <- sample(seq_along(seq_vec), n_mut)
  
  for (p in pos) {
    seq_vec[p] <- sample(bases[bases != seq_vec[p]], 1)
  }
  
  paste(seq_vec, collapse = "")
}

result_list <- list()

for(i in seq_along(top_sequence$sequence)){
  
  seq <- top_sequence$sequence[i]
  
  for(j in 1:10){
    
    n_mut <- sample(1:5,1)
    mut_seq <- mutate_seq(seq, n_mut)
    
    result_list[[length(result_list)+1]] <- data.frame(
      seq_id = i,
      mut_id = j,
      original = seq,
      mutated = mut_seq
    )
  }
}

mutated_df <- do.call(rbind, result_list)
head(mutated_df)
mutated_df$id <- paste0(mutated_df$seq_id,mutated_df$mut_id)
nrow(mutated_df)
write.csv(mutated_df,"mut.csv")
#################
mp <- read.csv("ensemble_probs_mut.csv")
ms <- read.csv("average_prob_importance_per_mut.csv")
mp <- mp[!duplicated(mp$sequence),]
top_sequence <- pos[1:10000,]
rownames(mp) <- mp$sequence
rownames(ms) <- ms$sequence
mp <- mp[rownames(ms),]
ms <- cbind(ms,mp)
ms <- ms[,-4]
ms$score <- ms$importance * ms$label_sum
ms$label <- "mutant"
top_sequence$label <- "raw"
colnames(ms)
colnames(top_sequence)
ms2 <- ms[,c(1,2,3,9,4,5,6,7,8,10)]
ms2 <- rbind(ms2,top_sequence)
ggplot(ms2, aes(x = score, fill = factor(label))) +
  geom_histogram(bins = 30, color = "black", alpha = 0.7) +
  labs(title = "Histogram of Score by Label", 
       x = "Score", 
       y = "Frequency") +
  scale_fill_manual(values = c("skyblue", "orange")) + 
  theme_minimal()
ms2 <- ms2[order(ms2$score,decreasing = T),]
write.csv(ms2[1:100,],"ProteinX.csv")