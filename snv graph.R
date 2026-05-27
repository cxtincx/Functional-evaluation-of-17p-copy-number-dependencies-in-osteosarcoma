
library(dplyr)
library(ggplot2)
library(stringr)
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)
library(rtracklayer)   
library(GenomicRanges)  
library(GenomeInfoDb)

snv_final_df$Start_Position <- as.numeric(as.character(snv_final_df$Start_Position))
snv_final_df$End_Position   <- as.numeric(as.character(snv_final_df$End_Position))

session <- browserSession("UCSC")
genome(session) <- "hg38"
query <- ucscTableQuery(session, table = "cytoBand")
cyto <- getTable(query)

cyto17p <- cyto %>%
  filter(chrom == "chr17", str_starts(name, "p")) %>%
  filter(gieStain != "acen")
cyto17p$chromStart <- as.numeric(cyto17p$chromStart)
cyto17p$chromEnd   <- as.numeric(cyto17p$chromEnd)

gie_cols <- c(gneg="white", gpos25="grey85", gpos50="grey70", gpos75="grey50", gpos100="grey30", acen="red")
cyto17p$fill_col <- gie_cols[ as.character(cyto17p$gieStain) ]

start_p <- min(cyto17p$chromStart, na.rm = TRUE)
end_p   <- max(cyto17p$chromEnd,  na.rm = TRUE)

snv_17p <- snv_final_df %>%
  filter(Start_Position >= start_p,
         End_Position   <= end_p,
         Chromosome %in% c("chr17"))

library(dplyr)
library(ggplot2)
have_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)


if (!exists("snv_17p")) stop("snv_17p not found.")
if (!exists("cyto17p")) stop("cyto17p not found.")
if (!exists("start_p") || !exists("end_p")) stop("start_p or end_p not found.")
if (!exists("Mart_chr17")) stop("Mart_chr17 not found.")

possible_names <- c("mut_class", "mutation_class", "Variant_Classification", "variant_class", "mutClass")
mut_col_name <- intersect(possible_names, colnames(snv_17p))
if (length(mut_col_name) == 0) {
  stop("No mutation-class column found in snv_17p among the usual candidates. Provide the column name.")
}
mut_col_name <- mut_col_name[1]  

mut_types_present <- unique(as.character(snv_17p[[mut_col_name]]))
cols_named <- setNames(
  grDevices::hcl(
    h = seq(15, 375, length = length(mut_types_present) + 1)[1:length(mut_types_present)],
    l = 65, c = 100
  ),
  mut_types_present
)


genes_df <- Mart_chr17 %>%
  filter(grepl("^17p", location)) %>%
  dplyr::select(gene) %>%
  distinct() %>%
  inner_join(
    snv_17p %>% group_by(Hugo_Symbol) %>% summarise(min_pos = min(Start_Position, na.rm = TRUE), .groups = "drop"),
    by = c("gene" = "Hugo_Symbol")
  ) %>%
  arrange(min_pos) %>%
  mutate(mid = as.numeric(min_pos), gene = as.character(gene)) %>%
  filter(mid >= start_p & mid <= end_p)

karyo_ymin  <- -0.25
karyo_ymax  <- 0.05
tick_ymax   <- -0.30  
tick_ymin   <- -0.65 
label_y_pos <- -0.70 
karyo_box <- data.frame(xmin = start_p, xmax = end_p, ymin = karyo_ymin, ymax = karyo_ymax)
snv_labels <- snv_17p %>%
  dplyr::group_by(Hugo_Symbol) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

p_final <- ggplot() +

  geom_rect(data = karyo_box,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = NA, color = "black", size = 0.45, inherit.aes = FALSE) +
  
  geom_rect(data = cyto17p,
            aes(xmin = chromStart, xmax = chromEnd, ymin = karyo_ymin, ymax = karyo_ymax),
            fill = cyto17p$fill_col, color = "black", size = 0.18, inherit.aes = FALSE) +
  geom_text(data = cyto17p,
            aes(x = (chromStart + chromEnd) / 2, y = (karyo_ymin + karyo_ymax) / 2, label = name),
            size = 4, fontface = "bold", inherit.aes = FALSE) +
 
  geom_segment(data = snv_17p,
               aes_string(x = "Start_Position", xend = "Start_Position",
                          y = as.character(tick_ymin), yend = as.character(tick_ymax),
                          color = mut_col_name),
               size = 0.7, linewidth = 1, inherit.aes = FALSE, show.legend = TRUE) +
 
 
  ggrepel::geom_text_repel(data = snv_labels,
                           aes(x = Start_Position, y = label_y_pos, label = Hugo_Symbol),
                           nudge_y = 0, angle = 90, segment.size = 0.25, segment.color = NA, 
                           hjust = 1, box.padding = 0.1, min.segment.length = 0.02, 
                           size = 4, inherit.aes = FALSE) +
  scale_color_manual(values = cols_named, name = "Mutation type") +
  coord_cartesian(xlim = c(start_p, end_p), ylim = c(-2, 1.2), expand = FALSE) +
  labs(x = "Chr17p SNVs", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 11),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.16, 0.1),
    legend.background = element_blank(),
    legend.key = element_rect(fill = NA),
    legend.title = element_text(size = 13, face = "bold"), 
    legend.text = element_text(size = 13),
    panel.background = element_rect(fill = "white", colour = NA),
    axis.text.x = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))



print(p_final)
ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 analysis plots /SNV_15052026.png", plot = p_final, width = 6, height = 4.5, dpi = 300)


