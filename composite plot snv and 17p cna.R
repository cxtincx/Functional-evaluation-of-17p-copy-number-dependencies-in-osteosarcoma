if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
library(ggplot2)
library(dplyr)

df <- OS1_copy_number_outrem[OS1_copy_number_outrem$chromosome == 17,]

df <- df %>% 
  filter(chromosome == "17" | chromosome == "chr17") %>% 
  mutate(chromosome = droplevels(factor(chromosome)))

df <- df %>%
  filter(chromosome %in% c("17")) %>% 
  arrange(order) %>%                    
  mutate(chr17_order = factor(row_number()))
         
average_data_OS1 <- df %>%
  group_by(chr17_order) %>%
  summarise(average_cn = mean(cn), .groups = 'drop')


positions_OS1_17p <- positions_OS1 %>%
  filter(chromosome %in% c("17", "chr17")) %>%
  mutate(
    largest_order_OS1 = nrow(df),
    mid = nrow(df) / 2
  )
df <- df %>%
  mutate(color = case_when(
    cn > 2 ~ "red",      # Red for cn > 2
    cn < 2 ~ "blue",     # Blue for cn < 2
    cn == 2 ~ "darkgrey"     # Grey for cn = 2
  )) 

df <- df %>% 
  left_join(chr17_OS1_copy_number_outrem %>% dplyr::select(gene, location), 
            by = c("gene" = "gene"))
df2 <- df %>%
  filter(location %in% c("17p12", "17p11.2", "17p13.3", "17p13.1", "17p11.2", "17p13", "17p13.2", "17p13.1-p12")) %>% 
  arrange(order) %>%                    
  mutate(chr17p_order = factor(row_number()))

df2 <- df2 %>% 
  left_join(OS1_CNV_calling_coding_filtered_unique_filtered %>% dplyr::select(coding_new, start, end), 
            by = c("gene" = "coding_new"))
df2 <- df2 %>% distinct(gene, .keep_all = TRUE)
average_data_OS1 <- average_data_OS1[base::order(as.numeric(average_data_OS1$chr17_order)), ]

rle_lengths <- rle(average_data_OS1$average_cn)$lengths
rle_values  <- rle(average_data_OS1$average_cn)$values
clean_segments <- data.frame(
  average_cn = rle_values,
  end_idx    = cumsum(rle_lengths)
)
clean_segments$start_idx <- c(1, clean_segments$end_idx[-nrow(clean_segments)] + 1)


OS1_cn_17<-ggplot(data = df2) +
  geom_point(aes(x=chr17_order, y = cn, color = color), size = 2, position = position_jitter(width = 0.5, height = 0.5),alpha = 1) +
 
  geom_segment(data = clean_segments, aes(x = chr17_order, xend = chr17_order, y = average_cn, yend = average_cn), color = "black", size = 0.5, linetype = "solid") +
  geom_text(data = positions_OS1, aes(x = mid, y = -0.7, label = chromosome), 
            vjust = 1, hjust = 0.5, color = "black", size = 7) + 
  geom_step(data = average_data_OS1, 
            aes(x = as.numeric(chr17_order), y = average_cn), 
            color = "black", 
            linetype = "solid", 
            linewidth = 0.8, 
            inherit.aes = FALSE) +
  theme_classic() + 
  labs(x = "Chromosome 17p", y = " Gene copy number", title = "NRH-OS1 copy number plot") +  
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size = 18))+
  theme(axis.title.x = element_text(size = 20))+
  theme(axis.title.y = element_text(size = 20))+
  theme(axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0.05, 0.05), breaks = seq(0, 25, by = 1))+
  scale_x_discrete(expand = c(0, 0.05))+  
  coord_cartesian(xlim = range(as.numeric(df2$chr17p_order), na.rm = TRUE)) +
  scale_color_identity()

print(OS1_cn_17)


ggsave("/Users/catincaapostol/Desktop/Master's/project/OS1_17p_cn_plot.png", OS1_cn_17, width = 7.5, height = 4.5, dpi = 300)

##karyo band 

library(dplyr)
library(rtracklayer)
library(stringr)
session <- browserSession("UCSC")
genome(session) <- "hg38"
query <- ucscTableQuery(session, table = "cytoBand")
cyto_raw <- getTable(query)
cyto17p <- cyto_raw %>%
  filter(chrom == "chr17", str_starts(name, "p")) %>%
  filter(gieStain != "acen") 
min_bp <- min(df2$start_position, na.rm = TRUE)
max_bp <- max(df2$end_position, na.rm = TRUE)

min_index <- 1
max_index <- max(as.numeric(df2$chr17p_order), na.rm = TRUE)

cyto17p_scaled <- cyto17p %>%
  mutate(
    scaled_start = min_index + (chromStart - min_bp) / (max_bp - min_bp) * (max_index - min_index),
    scaled_end   = min_index + (chromEnd - min_bp) / (max_bp - min_bp) * (max_index - min_index)
  )

gie_cols <- c(gneg="white", gpos25="grey85", gpos50="grey70", gpos75="grey50", gpos100="grey30", acen="red", gvar="grey90")
cyto17p_scaled$fill_col <- gie_cols[ as.character(cyto17p_scaled$gieStain) ]
cyto17p_scaled$fill_col[is.na(cyto17p_scaled$fill_col)] <- "white"
karyo_box <- data.frame(
  xmin = min_index,
  xmax = max_index,
  ymin = -1.5,  
  ymax = -0.5
)

average_data_OS1 <- average_data_OS1 %>%
  arrange(order) %>%
  mutate(chr17p_order = factor(row_number()))

OS1_cn_17 <- ggplot() +
  geom_point(data = df2, aes(x = chr17p_order, y = cn, color = color), 
             size = 2, position = position_jitter(width = 0.5, height = 0.5), alpha = 0.6) +
  
  geom_segment(data = average_data_OS1, aes(x = chr17p_order, xend = chr17p_order, y = average_cn, yend = average_cn), color = "black", size = 0.5) +
  geom_point(data = average_data_OS1, aes(x = chr17p_order, y = average_cn), color = "black", size = 0.001, shape = 0.01, fill = "black") + 
  geom_vline(data = positions_OS1, aes(xintercept = largest_order_OS1), linetype = "dashed", color = "grey", size = 0.5) + 
  
geom_rect(data = karyo_box,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
          fill = NA, color = "black", size = 0.5, inherit.aes = FALSE) +
  

  geom_rect(data = cyto17p_scaled,
            aes(xmin = scaled_start, xmax = scaled_end, ymin = -1.5, ymax = -0.5),
            fill = cyto17p_scaled$fill_col, color = "black", size = 0.2, inherit.aes = FALSE) +

  geom_text(data = cyto17p_scaled,
            aes(x = (scaled_start + scaled_end) / 2, y = -1.0, label = name),
            size = 3.5, fontface = "bold", inherit.aes = FALSE) +

theme_classic() +  
  labs(x = "Chromosome 17p", y = "Gene copy number", title = "NRH-OS1 copy number plot") +  
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(), 
    axis.text.y  = element_text(size = 14),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  ) +
  
  scale_y_continuous(breaks = seq(0, 25, by = 1), expand = c(0.02, 0.02)) +
  scale_x_discrete(expand = c(0, 0.05)) +  
  
  
  coord_cartesian(
    xlim = range(as.numeric(df2$chr17p_order), na.rm = TRUE),
    ylim = c(-2.0, 25)
  ) +
  scale_color_identity()

print(OS1_cn_17)

