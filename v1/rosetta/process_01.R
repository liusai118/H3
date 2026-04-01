data <- read.csv("combined_5mer_analysis_top80pct.csv")
data <- data[order(data$rank_mean,decreasing = F),]
kmers <- data$kmer
top_kmers <- kmers[1:100]
bot_kmers <- kmers[925:1024]
###################
library(ggplot2)
p1 <- ggplot(data, aes(
  x = rank_mean,
  y = log10(total_occurrences_sum),
  color = log10(total_occurrences_sum)
)) +
  geom_point(alpha = 0.5, size = 0.25,shape = 16) +
  scale_color_gradientn(
    colors = colorRampPalette(c(
      "#4575B4","#74ADD1","#ABD9E9","#E0F3F8",
      "#FFFFBF",
      "#FEE090","#FDAE61","#F46D43","#D73027"
    ))(100)
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Mean Rank",
    y = "log10(Total Frequency)",
    color = "log10(Frequency)",
    title = ""
  ) + 
  theme_classic()
p1
library(svglite)
ggsave(plot = p1, "p1.png",width =4 ,height = 2.5,dpi=2080)
###############################################

###############################################
library(data.table)
pos <- as.data.frame(fread("pos.csv"))
colnames(pos)
p2 <- ggplot(pos, aes(x = prob_mean, y = importance)) +
  geom_hex(bins = 200) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(
    x = "Mean Probability",
    y = "Importance",
    fill = "Count"
  ) + theme_classic()

ggsave(plot = p2, "p2.png",width =4 ,height = 2.5,dpi=2080)


########################################
library(data.table)
neg <- as.data.frame(fread("../sequence/neg_c.csv"))

p2 <- ggplot(neg, aes(x = prob_mean, y = importance)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(
    x = "Mean Probability",
    y = "Importance",
    fill = "Count"
  ) + theme_classic()
p2

p3 <- ggplot(pos, aes(x = prob_mean, y = importance)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c() +
  theme_bw() +
  labs(
    x = "Mean Probability",
    y = "Importance",
    fill = "Count"
  ) + theme_classic()
library(patchwork)
p2 | p3

ggsave(plot = p3, "p3.png",width =4 ,height = 2.5,dpi=2080)
###########################

