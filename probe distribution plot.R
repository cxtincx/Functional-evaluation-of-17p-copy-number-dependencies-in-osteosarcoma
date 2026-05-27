library(dplyr)  
library(ggplot2)        
library(readxl)
library(tidyr)
library(ggpubr)
library(writexl)

counts <- read_excel("/Users/catincaapostol/Desktop/Master's/project/os1 fish/fish counts.xlsx")
counts_long <- pivot_longer(counts,
                            cols = c("Green (TP53)", "Red (NCOR1)", "Aqua (centromere)"),
                            names_to = "group",
                            values_to = "values")

p <- ggplot(counts_long, aes(x = group, y = values, fill = group, colour = group)) +
  geom_violin(fill = NA, trim = FALSE, scale = "width", colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = c(
    "Green (TP53)" = "#4ff07a",
    "Red (NCOR1)" = "#f04f4f",
    "Aqua (centromere)" = "#4ff0eb"
  )) +
  geom_jitter(width = 0.1, height = 0, size = 4, shape = 16, alpha = 1) +
  stat_summary(fun = mean,
               fun.min = mean,
               fun.max = mean,
               geom = "crossbar",
               width = 0.5,
               fatten = 1,
               colour = "black",
               size = 0.7) +
  scale_colour_manual(values = c(
    "Green (TP53)" = "#4ff07a",
    "Red (NCOR1)" = "#f04f4f",
    "Aqua (centromere)" = "#4ff0eb"
  )) +
  stat_compare_means(comparisons = list(
    c("Green (TP53)", "Aqua (centromere)"), 
    c("Red (NCOR1)", "Aqua (centromere)"),
    c("Green (TP53)", "Red (NCOR1)")
  ),
  method = "wilcox.test",
  p.adjust.method = "bonferroni",
  aes(label = after_stat(p.signif)),
  size = 5,           
  bracket.size = 0.5,
  fontface = "bold") + 
  theme_classic() +
  labs(title = "", x = "", y = "Value", colour = "Group") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.75), 
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 18, colour = "black"), 
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

print (p)

stats_table <- compare_means(values ~ group, data = counts_long, 
                             method = "wilcox.test", 
                             p.adjust.method = "bonferroni")
print(stats_table)
write_xlsx(stats_table, "/Users/catincaapostol/Desktop/Master's/project/os1 fish/fish_stats_results.xlsx")
ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 fish/fish counts plot + statistical.pdf", plot = p, width = 6.5, height = 4.5, dpi = 300)


#chromosome counts 
chr_counts <- read_excel("/Users/catincaapostol/Desktop/Master's/project/os1 fish/chromosome counts .xlsx")
chr_counts_long <- pivot_longer(chr_counts,
                            cols = c("Chromosomes"),
                            names_to = "group",
                            values_to = "values")
p <- ggplot(chr_counts_long, aes(x = group, y = values, fill = group, colour = group, show.legend = FALSE)) +
  geom_violin(fill = NA, trim = FALSE, scale = "width", colour = "black") +
  geom_jitter(width = 0.1, height = 0, size = 4, shape = 16, alpha = 1) +
  stat_summary(fun = mean,
               fun.min = mean,
               fun.max = mean,
               geom = "crossbar",
               width = 0.5,
               fatten = 0,
               colour = "black",
               size = 0.7) +
  scale_colour_manual(values = c(
    "Chromosomes" = "darkblue"
  )) +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)))+
  labs(title = "", x = "", y = "Value", colour = "darkblue") +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.75), 
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 18, colour = "black"), 
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 18)
  )

print (p)
ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 fish/count chromosome plot + statistical.pdf", plot = p, width = 6.5, height = 4.5, dpi = 300)


#colocalization 


colocalized_counts <- read_excel("/Users/catincaapostol/Desktop/Master's/project/os1 fish/colocalization counts .xlsx")
colocalized_counts_long <- pivot_longer(colocalized_counts,
                            cols = c("Green + red", "Green + aqua", "Red + aqua"),
                            names_to = "group",
                            values_to = "values")

p <- ggplot(colocalized_counts_long, aes(x = group, y = values, fill = group, colour = group)) +
  geom_violin(fill = NA, trim = FALSE, scale = "width", colour = "black", show.legend = FALSE) +
  scale_colour_manual(values = c(
    "Green + red" = "darkgreen",
    "Green + aqua" = "darkblue",
    "Red + aqua" = "darkred"
  )) +
  geom_jitter(width = 0.1, height = 0, size = 4, shape = 16, alpha = 1) +
  stat_summary(fun = mean,
               fun.min = mean,
               fun.max = mean,
               geom = "crossbar",
               width = 0.5,
               fatten = 1,
               colour = "black",
               size = 0.7) +
  scale_colour_manual(values = c(
    "Green + red" = "darkgreen",
    "Green + aqua" = "darkblue",
    "Red + aqua" = "darkred"
  )) +
  stat_compare_means(comparisons = list(
    c("Green + red", "Green + aqua"), 
    c("Red + aqua", "Green + aqua"),
    c("Green + red", "Red + aqua")
  ),
  method = "wilcox.test",
  p.adjust.method = "bonferroni",
  aes(label = after_stat(p.signif)),
  size = 5,           
  bracket.size = 0.5,
  fontface = "bold") + 
  theme_classic() +
  labs(title = "", x = "", y = "Value", colour = "Group") +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.75), 
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(size = 18, colour = "black"), 
    axis.text.y = element_text(size = 18, colour = "black"),
    axis.title.y = element_text(size = 18),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )

print (p)

stats_table <- compare_means(values ~ group, data = colocalized_counts_long, 
                             method = "wilcox.test", 
                             p.adjust.method = "bonferroni")
print(stats_table)
write_xlsx(stats_table, "/Users/catincaapostol/Desktop/Master's/project/os1 fish/coloc_stats_results.xlsx")
ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 fish/colocalized plot.png", plot = p, width = 7, height = 4.5, dpi = 300)




