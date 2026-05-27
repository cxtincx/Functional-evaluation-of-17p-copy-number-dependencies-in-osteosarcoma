##log-fold changes for counts 
library(dplyr)
counts_df <- read.csv("/Users/catincaapostol/Desktop/Master's/project/chronos/OS1_count_matrix.csv")
library(dplyr)

library(dplyr)
counts_df <- counts_df[-27077,]
counts_df <- counts_df[-9429,]
normalized_df <- counts_df %>%
  mutate(
    T0_cpm = (T0 / sum(T0)) * 1e06,
    T1_cpm = (T1 / sum(T1)) * 1e06,
    T3_cpm = (T3 / sum(T3)) * 1e06,
    LFC = log2((counts_df$T3 + 1) / (counts_df$T0 + 1)))


negative_inzolia <- read.csv("/Users/catincaapostol/Desktop/Master's/project/chronos/Negative_controls_final2.csv")
negative_inzolia <- negative_inzolia %>%
  dplyr::left_join(
    dplyr::select(normalized_df, sgrna, LFC), 
    by = "sgrna"
  )
library(dplyr)


essential_genes <- read.csv("/Users/catincaapostol/Desktop/Master's/project/chronos/Hart_essential_genes_final.csv")

lfc_df <- normalized_df %>%
  mutate(sgrna_clean = sub("_[0-9]+$", "", sgrna)) %>%
  group_by(sgrna_clean) %>%
  summarise(LFC = mean(LFC), na.rm = TRUE)

essential_genes <- essential_genes %>%
  left_join(
    lfc_df, 
    by = c("gene" = "sgrna_clean")
  )
essential_genes <- essential_genes %>% 
  dplyr::select(-na.rm)

library(dplyr)

sgrna_guides <- normalized_df %>%
  dplyr::select(sgrna, LFC, T0_cpm, T3_cpm)


#plot

library(dplyr)
library(ggplot2)

##bottom 15% calculation
ranked_genes <- sgrna_guides %>%
  filter(is.finite(LFC)) %>%
  arrange(LFC)

total_genes <- nrow(ranked_genes)
cutoff_row <- round(total_genes * 0.15)

bottom_15 <- ranked_genes %>% 
  slice(1:cutoff_row)

bottom_15_genes <- sub("_[12]$", "", bottom_15$sgrna)
intersecting_genes <- intersect(bottom_15_genes, essential_genes$gene)
count_essentials <- length(intersecting_genes)
percent <- (count_essentials * 100 / (3152/2))
print(percent)


library(dplyr)

plot_df <- normalized_df %>%
  mutate(
    base_gene = sub("_\\d+$", "", sgrna),
    Group = case_when(
      sgrna %in% negative_inzolia$sgrna     ~ "Negative Controls",
      base_gene %in% essential_genes$gene   ~ "Essential Genes",
      TRUE                                  ~ "Library sgRNA"
    )
  ) %>%
  dplyr::select(-base_gene)


plot_df$Group <- factor(plot_df$Group, 
                          levels = c("Library sgRNA", "Negative Controls", "Essential Genes"))
p <- ggplot(plot_df, aes(x = LFC, fill = Group, color = Group)) +
  geom_density(alpha = 0.3, linewidth = 1) +
  scale_fill_manual(values = c("Library sgRNA" = "blue", 
                               "Negative Controls"   = "grey40", 
                               "Essential Genes"     = "red")) +
  scale_color_manual(values = c("Library sgRNA" = "blue", 
                                "Negative Controls"   = "grey40", 
                                "Essential Genes"     = "red")) +
  
  theme_classic(base_size = 20) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_text(size = 17), 
    legend.text = element_text(size = 17), 
    plot.caption = element_text(
      size = 16, 
      hjust = 0.5,       
      vjust = -1,   
      color = "black"
    ),
    plot.margin = margin(10, 10, 20, 10) 
  ) +
  labs(
    x = "Log2 Fold Change",
    y = "Density",
    fill = "Gene Category:",
    color = "Gene Category:", 
    caption = paste0("Essential Gene Recall: ", round(percent, 1), "%")
  ) 

print(p)
ggsave(filename = "/Users/catincaapostol/Desktop/Master's/project/chronos/lfc_plot.png", 
       plot = p, height = 5, width = 10, dpi = 300)

#z scores 
neg_control_lfcs <- plot_df %>% 
  filter(Group == "Negative Controls") %>% 
  pull(LFC)

mu_neg  <- mean(neg_control_lfcs, na.rm = TRUE)
sig_neg <- sd(neg_control_lfcs, na.rm = TRUE)

cassette_stats_df <- plot_df %>%
  mutate(
    cassette_z = (LFC - mu_neg) / sig_neg,
    cassette_p = pnorm(cassette_z),
    cassette_fdr = p.adjust(cassette_p, method = "BH")
  )

#keep genes in 17p12-p11.2

ordered_df <- amp_genes %>%
  left_join(
    cassette_stats_df %>% mutate(clean_name = str_remove(sgrna, "_[0-9]+$")),
    by = c("gene" = "clean_name")
  ) %>%
  filter(!is.na(gene))

#plot
library(ggplot2)
library(ggrepel)

ordered_df$gene <- factor(ordered_df$gene, 
                                 levels = unique(ordered_df$gene))

p <- ggplot(ordered_df, aes(x = gene, y = cassette_z)) +
  geom_point(color = "black", alpha = 1, size = 4) +
  theme_minimal(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20), 
    axis.text.y = element_text(color = "black")
  ) +
  labs(
    x = "Gene (Ordered)",
    y = "Z-Score"
  )
print(p)

ggsave(filename = "/Users/catincaapostol/Desktop/Master's/project/chronos/zscore.png", 
       plot     = p, height = 5, width = 11, dpi = 300)

