library(ggplot2)
library(Biostrings)
library(dplyr)
library(irlba)
library(Rtsne)
library(ggplot2)
library(hexbin)
library(Matrix)
all <- read.csv("../Dataset/rosetta.csv")
a <-  na.omit(all)
p7 <- ggplot(a, aes(log10(dSASA_int + 1), log10(abs(dG_separated) + 1))) +
  
  geom_point(
    data = subset(all, dG_separated <= 0),
    aes(color = fa_elec),
    alpha = 0.35,
    shape = 16,
    size = 2
  ) +
  
  geom_point(
    data = subset(all, dG_separated > 0),
    aes(color = fa_elec),
    shape = 4,
    alpha = 0.08,
    size = 1.8,
    stroke = 0.3
  ) +
  scale_color_viridis_c(
    option = "plasma",
    direction = -1,
    oob = scales::squish
  ) +
  
  guides(
    color = guide_colorbar(
      barheight = unit(5, "cm"),
      barwidth = unit(0.5, "cm")
    )
  ) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "log10 Interface area",
    y = "log10 Binding energy",
    color = "Electrostatic\nenergy"
  ) +
  
  theme(
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 11),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA)
  ) + 
  geom_rug(
    data = subset(all, dG_separated <= 0),
    color = "grey80",
    alpha = 0.15
  )
p7
ggsave(plot = p7,"t4.png",width = 5.5,height = 4.3,dpi = 2080,bg = "transparent")
####################################
b <- read.csv("../sequence/b.csv")
b$id <- b$id

sample_mean <- b %>%
  group_by(id, sample) %>%
  summarise(
    dG_separated = mean(dG_separated, na.rm = TRUE),
    dG_cross = mean(dG_cross, na.rm = TRUE),
    fa_elec = mean(fa_elec, na.rm = TRUE),
    delta_unsatHbonds = mean(delta_unsatHbonds, na.rm = TRUE),
    fa_intra_rep = mean(fa_intra_rep, na.rm = TRUE),
    dSASA_polar = mean(dSASA_polar, na.rm = TRUE),
    dSASA_int = mean(dSASA_int, na.rm = TRUE),
    hbond_E_fraction = mean(hbond_E_fraction, na.rm = TRUE),
    .groups = "drop"
  )

# 2 计算排名
sample_rank <- sample_mean %>%
  mutate(
    rank_dG_separated = rank(dG_separated, ties.method = "min"),
    rank_dG_cross = rank(dG_cross, ties.method = "min"),
    rank_elec = rank(fa_elec, ties.method = "min"),
    rank_unsat = rank(delta_unsatHbonds, ties.method = "min"),
    rank_rep = rank(fa_intra_rep, ties.method = "min"),
    
    rank_dSASA_polar = rank(-dSASA_polar, ties.method = "min"),
    rank_dSASA_int = rank(-dSASA_int, ties.method = "min"),
    rank_hbond = rank(-hbond_E_fraction, ties.method = "min")
  )

# 3 计算综合排名
sample_rank <- sample_rank %>%
  mutate(
    total_rank = rowMeans(
      select(., starts_with("rank_"))
    )
  )

# 4 每个 complex 选最佳构象
best_complex <- sample_rank %>%
  group_by(id) %>%
  slice_min(order_by = total_rank, n = 1, with_ties = FALSE) %>%
  ungroup()

best_complex <- best_complex %>%
  arrange(total_rank)

b_rank <- best_complex %>%
  group_by(id) %>%
  summarise(
    dG_separated = mean(dG_separated, na.rm = TRUE),
    dG_cross = mean(dG_cross, na.rm = TRUE),
    fa_elec = mean(fa_elec, na.rm = TRUE),
    delta_unsatHbonds = mean(delta_unsatHbonds, na.rm = TRUE),
    fa_intra_rep = mean(fa_intra_rep, na.rm = TRUE),
    dSASA_int = mean(dSASA_int, na.rm = TRUE),
    hbond_E_fraction = mean(hbond_E_fraction, na.rm = TRUE)
  ) %>%
  mutate(
    rank_dG_separated = rank(dG_separated, ties.method = "min"),
    rank_dG_cross = rank(dG_cross, ties.method = "min"),
    rank_elec = rank(fa_elec, ties.method = "min"),
    rank_unsat = rank(delta_unsatHbonds, ties.method = "min"),
    rank_rep = rank(fa_intra_rep, ties.method = "min"),
    rank_dSASA_int = rank(-dSASA_int, ties.method = "min"),
    rank_hbond = rank(-hbond_E_fraction, ties.method = "min")
  )
b_rank
b_rank <- b_rank %>%
  mutate(
    total_rank = rowMeans(select(., starts_with("rank_")))
  ) %>%
  arrange(total_rank)

b_rank <- read.csv("rank1.csv")
b_rank <- b_rank[order(b_rank$total_rank),]
df10 <- b_rank[1:10,]

write.csv(df10,"result.csv")