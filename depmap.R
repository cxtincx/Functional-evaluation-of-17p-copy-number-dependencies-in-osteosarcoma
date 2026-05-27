library(dplyr)
library(ggplot2)
library(stringr)
library(ggrepel)
library(rtracklayer)   
library(GenomicRanges)  
library(GenomeInfoDb)

crispr_gene_effect <- read.csv("/Users/catincaapostol/Desktop/Master's/project/depmap/DepMap_dfs/CRISPRGeneEffectUncorrected.csv")
screen_gene_effect <- read.csv("/Users/catincaapostol/Desktop/Master's/project/depmap/DepMap_dfs/ScreenGeneEffectUncorrected.csv")

screen_map <- read.csv("/Users/catincaapostol/Desktop/Master's/project/depmap/DepMap_dfs/CRISPRScreenMap.csv")
cn_screens <- read.csv("/Users/catincaapostol/Desktop/Master's/project/depmap/DepMap_dfs/Omics_Absolute_CN_Gene_Public_24Q4_subsetted.csv")
model <- read.csv("/Users/catincaapostol/Desktop/Master's/project/depmap/DepMap_dfs/Model.csv")
library_data <- read.csv("/Users/catincaapostol/Desktop/Master's/project/depmap/DepMap_dfs/ScreenSequenceMap.csv")

print(colnames(screen_map))
screen_gene_effect <- screen_gene_effect %>% 
  dplyr::rename(ScreenID = 1)

crispr_gene_effect <- crispr_gene_effect %>% 
  dplyr::rename(ModelID = 1)
cn_screens <- cn_screens %>% 
  dplyr::rename(ModelID = 1)
combined_df <- screen_gene_effect %>%
  left_join(screen_map, by = "ScreenID") %>%
  left_join(crispr_gene_effect, by = "ModelID", suffix = c("_screen", "_crispr"))

#filter for cn
crispr_gene_effect <- crispr_gene_effect %>%
  dplyr::filter(crispr_gene_effect$ModelID %in% cn_screens$ModelID)

screen_map_with_cn <- screen_map %>%
  dplyr::filter(screen_map$ModelID %in% crispr_gene_effect$ModelID)

screen_gene_effect <- screen_gene_effect %>%
  dplyr::filter(screen_gene_effect$ScreenID %in% screen_map_with_cn$ScreenID)

#filter for humagne 
library_data <- library_data[library_data$Library %in% c("Humagne-CD"),]
screen_map_final <- screen_map_with_cn %>%
  dplyr::filter(screen_map_with_cn$ModelID %in% library_data$ModelID)
cn_screens <- cn_screens %>%
  dplyr::filter(cn_screens$ModelID %in% screen_map_final$ModelID)
crispr_gene_effect <- crispr_gene_effect %>%
  dplyr::filter(crispr_gene_effect$ModelID %in% cn_screens$ModelID)
screen_map_cn_final <- screen_map_final %>%
  dplyr::filter(screen_map_final$ModelID %in% cn_screens$ModelID)
screen_gene_effect <- screen_gene_effect %>%
  dplyr::filter(screen_gene_effect$ScreenID %in% screen_map_cn_final$ScreenID)

#zoom in to amp genes
interest_genes <- amp_genes$gene
metadata_cols <- c("ModelID") 

combined_cn <- crispr_gene_effect %>%
  left_join(cn_screens, by = "ModelID", suffix = c("_crispr", "_cn"))
clean_column_names <- sub("\\.\\..*", "", colnames(combined_cn))
target_genes <- which(clean_column_names %in% interest_genes)

depmap_gene_columns <- colnames(combined_cn)[target_genes]
filtered_crispr_df <- combined_cn %>%
  dplyr::select(ModelID, all_of(depmap_gene_columns))

#find percentage of target genes amplified 

existing_genes <- intersect(interest_genes, colnames(filtered_crispr_df))
percentage_genes_amp <- filtered_crispr_df %>%
  rowwise() %>%
  mutate(
    amp_count = sum(c_across(any_of(existing_genes)) >= 3, na.rm = TRUE),
    percentage = (amp_count/length(existing_genes)) * 100
  ) %>%
  ungroup() %>%
  dplyr::select(ModelID, percentage)
head(percentage_genes_amp)
library(dplyr)

percentage_genes_amp <- percentage_genes_amp %>%
  mutate(
    group = case_when(
      percentage < 30  ~ "Non-amplified",
      percentage > 70  ~ "Amplified",
      TRUE             ~ "Low-level amplification"
    )
  )
filtered_crispr_df <- filtered_crispr_df %>%
  left_join(percentage_genes_amp, by = "ModelID")

#heatmap 
library(pheatmap)
library(dplyr)

ordered_genes_list <- amp_genes %>%
  arrange(`17p order`) %>% 
  pull(gene)

heatmap_data_prep <- filtered_crispr_df %>%
  filter(group %in% c("Amplified", "Non-amplified", "Low-level amplification")) %>%
  mutate(group = ifelse(group == "Low-level amplification", "Non-amplified", group)) %>%
  arrange(desc(group))

num_amplified <- sum(heatmap_data_prep$group == "Amplified")

all_chronos_cols <- colnames(heatmap_data_prep)[grepl("\\.\\.", colnames(heatmap_data_prep))]

chronos_gene_cols <- all_chronos_cols[sub("\\.\\..*", "", all_chronos_cols) %in% ordered_genes_list]

gene_name_prefixes <- sub("\\.\\..*", "", chronos_gene_cols)
chronos_gene_cols  <- chronos_gene_cols[match(intersect(ordered_genes_list, gene_name_prefixes), gene_name_prefixes)]

p_values <- sapply(chronos_gene_cols, function(gene_col) {
  amp_scores <- heatmap_data_prep[[gene_col]][heatmap_data_prep$group == "Amplified"]
  non_amp_scores <- heatmap_data_prep[[gene_col]][heatmap_data_prep$group == "Non-amplified"]
  amp_scores <- amp_scores[!is.na(amp_scores)]
  non_amp_scores <- non_amp_scores[!is.na(non_amp_scores)]
  if(length(amp_scores) > 1 & length(non_amp_scores) > 1) {
    return(wilcox.test(amp_scores, non_amp_scores)$p.value)
  } else {
    return(NA)
  }
})
p_adjusted <- p.adjust(p_values, method = "BH")

heatmap_matrix <- as.matrix(heatmap_data_prep[, chronos_gene_cols])
mode(heatmap_matrix) <- "numeric"

clean_names <- sub("\\.\\..*", "", colnames(heatmap_matrix))
colnames(heatmap_matrix) <- ifelse(!is.na(p_adjusted) & p_adjusted < 0.05, paste0(clean_names, "*"), clean_names)

row_labels <- rep("", nrow(heatmap_matrix))
row_labels[round(num_amplified / 2)] <- "AMPLIFIED MODELS"
row_labels[num_amplified + round((nrow(heatmap_matrix) - num_amplified) / 2)] <- "NON-AMPLIFIED MODELS"
rownames(heatmap_matrix) <- row_labels

palette_length <- 50
my_breaks <- c(
  seq(-1.5, -0.01, length.out = palette_length / 2), 
  0,                                                               
  seq(0.01, 0.5, length.out = palette_length / 2)  
)

p <- pheatmap(
  heatmap_matrix,
  cluster_rows = FALSE,            
  cluster_cols = FALSE,            
  show_rownames = FALSE,             
  show_colnames = TRUE,            
  gaps_row = num_amplified,        
  breaks = my_breaks,              
  na_col = "grey60",            
  color = colorRampPalette(c("blue", "beige", "red"))(palette_length),
  legend_labels = c("Highly Dependent", "No Effect (Neutral)", "Resistant"),
  fontsize_col = 17,
  fontsize = 16, 
  filename = "/Users/catincaapostol/Desktop/Master's/project/depmap/Amplification_Heatmap_stats.png",
  width = 15, 
  height = 6, 
  res = 300
)
print(p)
 

png("/Users/catincaapostol/Desktop/Master's/project/depmap/Amplification_Heatmap.png", width = 15, height = 4.5, dpi = 300)
stats_table <- data.frame(
  Gene = sub("\\.\\..*", "", chronos_gene_cols),
  Raw_P_Value = p_values,
  FDR_Corrected_P_Value = p_adjusted,
  row.names = NULL
  ) %>% 
  arrange(FDR_Corrected_P_Value)

write.csv(stats_table, "/Users/catincaapostol/Desktop/Master's/project/depmap/Gene_Dependency_Stats.csv", row.names = FALSE)

dev.off()

## p values

library(ComplexHeatmap)
library(dplyr)

ordered_matrix <- heatmap_matrix[, stats_table$Gene, drop = FALSE]
fdr_strings <- sprintf("%.2g", stats_table$FDR_Corrected_P_Value)

top_p_vals = HeatmapAnnotation(
  FDR = anno_text(fdr_strings, gp = gpar(fontsize = 14))
)

png(
  filename = "/Users/catincaapostol/Desktop/Master's/project/depmap/Amplification_Heatmap_stats.png", 
  width = 15, 
  height = 6, 
  units = "in", 
  res = 300
)

pheatmap(
  ordered_matrix,
  cluster_rows = FALSE,            
  cluster_cols = FALSE,            
  show_rownames = FALSE,             
  show_colnames = TRUE,            
  gaps_row = num_amplified,        
  breaks = my_breaks,              
  na_col = "grey60",            
  color = colorRampPalette(c("red", "beige", "blue"))(palette_length),
  
  heatmap_legend_param = list(title = "Chronos Score", title_gp = gpar(fontsize = 16, fontface = "bold"), 
                              labels_gp = gpar(fontsize = 16)  ),
  
  fontsize_col = 17,
  fontsize = 16,
  top_annotation = top_p_vals
)

dev.off()

#violin 
lowest_p <- stats_table$Gene[which.min(stats_table$FDR_Corrected_P_Value)]

PMP22 <- filtered_crispr_df %>% 
  select(starts_with(lowest_p))
colnames(PMP22) <- c("Chronos Score", "Copy Number")
PMP22 <- PMP22 %>%
  mutate(
    Amplification = if_else((`Copy Number`) > 2, "Amplified", "Non-amplified")
  )
plot_df <- PMP22 
library(ggplot2)
library(ggpubr) 

p <- ggplot(plot_df, aes(x = Amplification, y = `Chronos Score`, fill = Amplification)) +
  geom_violin(trim = FALSE, alpha = 0, show.legend = FALSE, color = "black") +
  geom_jitter(width = 0.15, shape = 16, alpha = 1, size = 4, aes(color = Amplification), show.legend = TRUE) +
  
  theme_classic(base_size = 20) +
  labs(
    x = "Amplification Status",
    y = "Chronos Score"
  ) +
  scale_fill_manual(values = c("Amplified" = "red", "Non-amplified" = "blue")) +
  scale_color_manual(values = c("Amplified" = "red", "Non-amplified" = "blue")) +
  stat_summary(fun = mean,
               fun.min = mean,
               fun.max = mean,
               geom = "crossbar",
               width = 0.5,
               fatten = 1,
               colour = "black",
               size = 0.7) +
  stat_compare_means (
    method = "wilcox.test", 
    label = "p.format",   
    label.x = 1.5,         
    size = 5,           
    bracket.size = 0.5,
    fontface = "bold") 

print(p)


 ggsave(
   filename = "/Users/catincaapostol/Desktop/Master's/project/depmap/PMP22_violin_TRUE.png", 
   plot = p,
   device = ragg::agg_png, 
   width = 7, 
   height = 5, 
   dpi = 300
 )




