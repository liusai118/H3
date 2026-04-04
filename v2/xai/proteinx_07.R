pos <- read.csv("../Dataset/pos.csv")
pos <- pos[order(pos$score,decreasing = T),]
top_sequence <- pos[1:10000,]
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