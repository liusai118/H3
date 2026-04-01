library(Biostrings)
library(dplyr)
library(irlba)
library(Rtsne)
library(ggplot2)
library(hexbin)
library(Matrix)
pos <- read.csv("pos.csv")
stopifnot("sequence" %in% colnames(pos))

set.seed(123)

sample_n <- min(400000, nrow(pos))

pos_sample <- pos %>%
  mutate(
    sequence = toupper(sequence),
    sequence = gsub("\\s+", "", sequence)
  ) %>%
  filter(
    grepl("^[ACGT]+$", sequence),
    nchar(sequence) >= 15,
    nchar(sequence) <= 20
  ) %>%
  slice_sample(n = sample_n)

dna_sample <- DNAStringSet(pos_sample$sequence)
k <- 5
kmer_mat <- oligonucleotideFrequency(
  dna_sample,
  width = k,
  step = 1,
  as.prob = TRUE,
  with.labels = TRUE
)
kmer_mat <- as.matrix(kmer_mat)
kmer_scaled <- scale(kmer_mat)
n_pcs <- 30
library(irlba)
library(RhpcBLASctl)
set.seed(123)
blas_set_num_threads(16)
omp_set_num_threads(16)


pca_res <- prcomp_irlba(
  kmer_scaled,
  n = n_pcs,
  center = FALSE,
  scale. = FALSE
)

pc_df <- as.data.frame(pca_res$x)
colnames(pc_df) <- paste0("PC", seq_len(ncol(pc_df)))

pca_var <- pca_res$sdev^2
pca_var_ratio <- pca_var / sum(pca_var)

print(round(pca_var_ratio[1:min(10, length(pca_var_ratio))], 4))


set.seed(123)
tsne_res <- Rtsne(
  as.matrix(pc_df),
  dims = 2,
  perplexity = 30,
  theta = 0.5,
  verbose = TRUE,
  max_iter = 1000,
  pca = FALSE,
  check_duplicates = FALSE,num_threads = 4
)

tsne_df <- data.frame(
  tSNE1 = tsne_res$Y[, 1],
  tSNE2 = tsne_res$Y[, 2]
)


plot_tsne_df <- bind_cols(pos_sample, pc_df, tsne_df)


library(viridis)

lens <- sort(unique(plot_tsne_df$true_len))

cols <- viridis(length(lens), option = "plasma")
names(cols) <- lens

p <- ggplot(plot_tsne_df, aes(tSNE1, tSNE2, color = factor(true_len))) +
  geom_point(size = 0.25, alpha = 0.6,shape = 16) +
  scale_color_manual(values = cols) +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Sequence length"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) + 
  theme(
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA)
  )

p
ggsave(plot = p,"t2.png",width = 6.5,height = 4.3,dpi = 2080,bg = "transparent")
#######################
pal_simpsons(palette = c("springfield"), alpha = 1)(10)

library(ggplot2)
ms2 <- rbind(neg,pos)
ggplot(ms2, aes(x = score, fill = factor(label))) +
  geom_histogram(
    bins = 35,
    color = "white",
    alpha = 0.80,
    position = "identity"
  ) +
  scale_fill_manual(
    values = c("#FED439FF", "#709AE1FF"),
    name = "Label"
  ) +
  labs(
    title = "",
    subtitle = "",
    x = "Score",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 14) +
  theme_classic(

  )
  ############################
  library(ggplot2)
#########################
pos <- read.csv("pos.csv")
l15 <- pos[pos$true_len == 15,]
l16 <- pos[pos$true_len == 16,]
l17 <- pos[pos$true_len == 17,]
l18 <- pos[pos$true_len == 18,]
l19 <- pos[pos$true_len == 19,]
l20 <- pos[pos$true_len == 20,]
mean(l15$score)
mean(l16$score)
mean(l17$score)
mean(l18$score)
mean(l19$score)
mean(l20$score)
##########################
table(pos$true_len)
ggplot(pos, aes(x = score)) +
  geom_histogram(bins = 25, fill = "#4DBBD5", color = "black") +
  theme_classic() +
  labs(x = "score", y = "Count")
library(ggplot2)

ggplot(pos, aes(x = score)) +
  geom_histogram(
    bins = 25,
    fill = "#4DBBD5",
    color = "white",
    alpha = 0.9
  ) +
  labs(
    title = "Distribution of Score",
    x = "Score",
    y = "Count"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major.y = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  )

ggplot(pos, aes(x = score)) +
  geom_histogram(aes(y = ..density..),
                 bins = 25,
                 fill = "#4DBBD5",
                 color = "white",
                 alpha = 0.8) +
  geom_density(color = "#D55E00", size = 1.2) +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of Score",
    x = "Score",
    y = "Density"
  )

##########################
library(dplyr)
pos %>%
  group_by(true_len) %>%
  summarise(
    mean_score = mean(score, na.rm = TRUE),
    sd = sd(score, na.rm = TRUE),
    n = n()
  )

###########################################
library(ggpubr)
ggplot(pos, aes(x = factor(true_len), y = score, fill = factor(true_len))) +
  geom_violin(trim = FALSE, alpha = 0.8, color = "grey30", linewidth = 0.3) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, linewidth = 0.3) +
  scale_fill_manual(values = c(
    "#C6DBEF", "#9ECAE1", "#6BAED6",
    "#4292C6", "#2171B5", "#084594"
  )) +
  theme_classic(base_size = 14) +
  labs(
    x = "True length",
    y = "Score",
    title = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  )
#########