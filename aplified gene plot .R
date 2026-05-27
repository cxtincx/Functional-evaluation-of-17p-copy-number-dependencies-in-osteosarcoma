#genes in region of amplification 

amplification <- chr17_OS1_copy_number_outrem[chr17_OS1_copy_number_outrem$location %in% c("17p12", "17p11.2"),]


amplification <- amplification[amplification$order_17 > 251,]

amp_genes <- amplification %>%  
  left_join(OS1_final_cn %>% 
              dplyr::select(start, end, gene),
            by = "gene"
            )
library(dplyr)
library(ggplot2)
library(stringr)
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
library(rtracklayer)   
library(GenomicRanges)  
library(GenomeInfoDb)

amp_genes$start <- as.numeric(as.character(amp_genes$start))
amp_genes$end   <- as.numeric(as.character(amp_genes$end))

amp_genes <- amp_genes %>%
  mutate(`17p order` = row_number())

band_data <- cyto17p_subset %>%
  mutate(
    band_label = paste0("17", name) 
  )

left_coord  <- min(amp_genes$start, na.rm = TRUE)
right_coord <- max(amp_genes$end, na.rm = TRUE)

library(ggplot2)
library(dplyr)

band_data <- cyto17p_subset %>%
  mutate(band_label = paste0("17", name))

left_coord  <- min(amp_genes$start, na.rm = TRUE)
right_coord <- max(amp_genes$end, na.rm = TRUE)

# 3. Pull the minimum and maximum order numbers to anchor the top coordinates
min_order <- min(amp_genes$`17p order`, na.rm = TRUE)
max_order <- max(amp_genes$`17p order`, na.rm = TRUE)

p <- ggplot() +
  geom_rect(data = band_data,
            aes(xmin = min_order - 0.5, xmax = max_order + 0.5, ymin = -0.2, ymax = 0.2, fill = fill_col),
            color = "black", linewidth = 0.5) +
  geom_text(aes(x = (min_order + max_order) / 2, y = 0, label = "17p12-17p11.2 amplified genes"),
            fontface = "bold", size = 7) +
  
  geom_text(aes(x = min_order, y = 0.5, label = as.character(left_coord)),
            hjust = 0, fontface = "bold", size = 5) +
  geom_text(aes(x = max_order, y = 0.5, label = as.character(right_coord)),
            hjust = 1, fontface = "bold", size = 5) +
  geom_segment(data = amp_genes,
               aes(x = `17p order`, xend = `17p order`, y = -0.2, yend = -0.5),
               linewidth = 0.5, color = "black") +
  geom_text(data = amp_genes,
            aes(x = `17p order`, y = -0.6, label = gene),
            angle = 90, hjust = 1, size = 4.5) +
  
  scale_fill_identity() +
  xlim(min_order - 1, max_order + 1) +
  ylim(-1.5, 1) +
  theme_void()

p <- p +
theme(
  plot.background = element_rect(fill = "white", color = NA)
)
print(p)


ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 analysis plots /amp_genes.png", plot = p, width = 14, height = 4.5, dpi = 300)





