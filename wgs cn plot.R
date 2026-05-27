
if (!requireNamespace("ggplot2", quietly=TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly=TRUE)) install.packages("dplyr")
library(ggplot2)
library(dplyr)

 df <- read.csv("/Users/catincaapostol/Desktop/Master's/project/wgs os1/OS1_copy_number.csv")

 df <- df %>%
   mutate(
     order = as.numeric(order),
     cn = as.numeric(cn),
     gene = as.character(gene),
     chromosome = as.character(chromosome)
   ) %>%
   arrange(order)
 
 # ---- categorize CN ----
 tol <- 0.1
 df2 <- df %>%
   mutate(cn_cat = case_when(
     is.na(cn)        ~ "missing",
     cn < 2 - tol     ~ "loss",
     cn > 2 + tol     ~ "gain",
     TRUE             ~ "neutral"
   ))
 
 # ---- chromosome boundaries ----
 chr_bounds <- df2 %>%
   group_by(chromosome) %>%
   summarize(
     x_min = min(order, na.rm=TRUE),
     x_max = max(order, na.rm=TRUE),
     x_mid = (min(order,na.rm=TRUE) + max(order,na.rm=TRUE)) / 2
   ) %>%
   arrange(x_min)
 
 # ---- per-chromosome median CN ----
 chr_median <- df2 %>%
   group_by(chromosome) %>%
   summarize(
     median_cn = median(cn, na.rm=TRUE),
     x_min = min(order, na.rm=TRUE),
     x_max = max(order, na.rm=TRUE)
   ) %>%
   arrange(x_min)
 
 cols <- c(loss="blue", neutral="gray50", gain="red", missing="black")
 
 # ---- PLOT ----
 p <- ggplot(df2, aes(x = order, y = cn)) +
   geom_point(aes(color = cn_cat), size = 0.8, alpha = 0.75) +
   geom_segment(data = chr_median,
                aes(x = x_min, xend = x_max, y = median_cn, yend = median_cn),
                inherit.aes = FALSE,
                color = "black",
                size = 0.9) +
   geom_vline(
     xintercept = if (nrow(chr_bounds) > 1) chr_bounds$x_min[-1] else NULL,
     linetype = "dashed",
     color = "grey80"
   ) +
   geom_hline(yintercept = 2, linetype = "dashed") +
   scale_x_continuous(
     breaks = chr_bounds$x_mid,
     labels = chr_bounds$chromosome,
     expand = c(0.01, 0.01)
   ) +
   scale_color_manual(values = cols, name = NULL) +
   labs(
     x = "Chromosome (genome-wide order)",
     y = "Copy Number",
     title = "Genome-wide CN Plot"
   ) +
   theme_minimal() +
   theme(axis.text.x = element_text(size = 10, face = "bold"))
 
 print(p)

# optionally save:
ggsave("OS1_wgs_cn_plot.png", p, width = 14, height = 4, dpi = 300)
