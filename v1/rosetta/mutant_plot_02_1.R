library(ggplot2)
library(Biostrings)
library(dplyr)
library(irlba)
library(Rtsne)
library(ggplot2)
library(hexbin)
library(Matrix)
## ================================
## ================================
all <- read.csv("../Dataset/mut.csv")
mut <- all[all$label == "mutant",]
dna_sample <- DNAStringSet(mut$sequence)

k <- 5

kmer_mat <- oligonucleotideFrequency(
  dna_sample,
  width = k,
  step = 1,
  as.prob = TRUE,
  with.labels = TRUE
)

kmer_mat <- as.matrix(kmer_mat)

n_pcs <- 30


library(irlba)
library(RhpcBLASctl)
set.seed(123)
blas_set_num_threads(32)
omp_set_num_threads(32)

pca_res <- prcomp_irlba(
  kmer_mat,
  n = 20,
  center = FALSE,
  scale. = FALSE,
  maxit = 1000,
  work = 80
)
pc_df <- as.data.frame(pca_res$x)
colnames(pc_df) <- paste0("PC", seq_len(ncol(pc_df)))


pca_var <- pca_res$sdev^2
pca_var_ratio <- pca_var / sum(pca_var)


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

plot_tsne_df <- bind_cols(mut, pc_df, tsne_df)

p <- ggplot(plot_tsne_df, aes(tSNE1, tSNE2, color = score)) +
  geom_point(size = 0.8, alpha = 1, shape = 16) +
  scale_color_gradientn(
    colors = colorRampPalette(c(
      "#4575B4","#74ADD1","#ABD9E9","#E0F3F8",
      "#FFFFBF",
      "#FEE090","#FDAE61","#F46D43","#D73027"
    ))(100)
  ) +
  labs(
    x = "t-SNE 1",
    y = "t-SNE 2",
    color = "Score"
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
ggsave(plot = p,"t3.png",width = 5.5,height = 4.3,dpi = 2080,bg = "transparent")
