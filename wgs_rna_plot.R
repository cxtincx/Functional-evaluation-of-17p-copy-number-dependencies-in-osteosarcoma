library(dplyr)

wgs_rna_df_common <- genes_keep_df %>%
  dplyr::select(
    symbol,
    padj,
    negLog10Padj,
    group,
    log2FoldChange,
    `Karyotype band.x`,
    `Chromosome/scaffold name.x`
  ) %>%
  dplyr::inner_join(
    OS1_copy_number %>%
      dplyr::select(gene, cn, chromosome),
    by = c("symbol" = "gene")
  ) %>%
  dplyr::filter(`Chromosome/scaffold name.x` == chromosome)
#filter for region of interest 
wgs_rna_df_common <- wgs_rna_df_common %>%
  dplyr::filter(
    symbol %in% amp_genes$gene
  )
#filter group to drop NS
wgs_rna_df_common <- wgs_rna_df_common %>%
  dplyr::filter(group != "NS") 
padj_cutoff <- 0.05

wgs_rna_df_common$group <- "NS"
wgs_rna_df_common$group[ wgs_rna_df_common$padj < padj_cutoff & wgs_rna_df_common$log2FoldChange > 0 ] <- "NRH-OS1"
wgs_rna_df_common$group[ wgs_rna_df_common$padj < padj_cutoff & wgs_rna_df_common$log2FoldChange < 0 ] <- "HOS"


library(ggplot2)
library(dplyr)

library(ggplot2)
library(scales) 
p <- ggplot(
  wgs_rna_df_common,
  aes(x = cn, y = log2FoldChange, color = group)
) +
  geom_point(size = 4, alpha = 1) +
  
  geom_text_repel(
    aes(label = symbol),
    size = 5,
    color = "black",
    max.overlaps = 5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("NRH-OS1" = "firebrick", "HOS" = "steelblue"),
    name = "Group"
  ) +

  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Gene copy number",
    breaks = pretty_breaks(n = 6), 
    limits = c(0,18)
  ) +
  scale_y_continuous(
    name = "log2Fold Change",
    limits = c(0, 2)
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = "NRH-OS1",
    hjust = 1.02, vjust = 1.2,
    color = "firebrick",
    size = 6
  ) +
  annotate(
    "text",
    x = Inf, y = -Inf,
    label = "HOS",
    hjust = 1.02, vjust = -0.3,
    color = "steelblue",
    size = 6
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black"),
    panel.grid.minor = element_blank(), 
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 20), 
    axis.text.x = element_text(size = 18), 
    axis.text.y = element_text(size = 18)
  )

print(p)
ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 analysis plots /17p_OS1andOHSvsHOS_CNvsRNA.png", p, width = 7, height = 4.5, dpi = 300)


