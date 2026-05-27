#mrna upreg genes in OHS and OS1 
#x axis - log2fold change 
#y axis - CN

library(ggplot2)
library(dplyr)
library(ggrepel)

df_both_lines <- rna_genes_both %>%
  filter(rna_genes_both$`Chromosome/scaffold name.x` == 17 & rna_genes_both$'Karyotype band.x' %in% c("p12", "p11.2")) %>%
  left_join(OS1_final_cn_order, by = c("symbol" = "gene"))%>%
  drop_na()

p <- ggplot(df_both_lines, aes(x = `cn`, y = `log2FoldChange`, label = `symbol`)) +
  geom_point(color = "firebrick", alpha = 0.6) + 
  geom_text_repel(
    nudge_x = 0.2,      
    direction = "y",
    hjust = 0,
    size = 3,
    color = "black"
  ) +
  scale_y_log10(labels = scales::comma) +
  theme_bw() +
  labs(
    title = "Genes Overexpressed in NRH-OS1 and OHS compared to HOS",
    x = "NRH-OS1 gene copy number",
    y = "log2Fold Change"
  )

 print(p)
 
ggsave("/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/OS1 and OHS overexpressed.png", width = 8, height = 3, dpi = 300)
 
 