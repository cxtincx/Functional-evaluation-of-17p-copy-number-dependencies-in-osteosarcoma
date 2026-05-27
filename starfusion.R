library(readr)
library(dplyr)
star_1_df <- read_tsv("/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/star fusion/1_finspector.FusionInspector.fusions.abridged.tsv")
star_2_df <- read_tsv("/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/star fusion/2_finspector.FusionInspector.fusions.abridged.tsv")
star_3_df <- read_tsv("/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/star fusion/3_finspector.FusionInspector.fusions.abridged.tsv")

combined <- bind_rows(
  star_1_df %>% mutate(Replicate = "R1"),
  star_2_df %>% mutate(Replicate = "R2"),
  star_3_df %>% mutate(Replicate = "R3")
)

names(combined)[1] <- "FusionName"
final_fusions <- combined %>%
  group_by(FusionName, LeftGene, RightGene, LeftBreakpoint, RightBreakpoint) %>%
  summarize(
    TimesDetected = n(),              
    MeanJunctionReads = mean(JunctionReadCount),
    MeanSpanningFrags = mean(SpanningFragCount),
    MeanFFPM = mean(FFPM),
    .groups = 'drop'
  ) %>%
  filter(TimesDetected >= 2, MeanJunctionReads >3, MeanFFPM >= 0.1) 
write_tsv(final_fusions, "/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/star fusion/final_fusions.tsv")

#plot 

png("/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/star fusion/Fusions NRH-OS1.png", 
    width = 5, height = 5, units = "in", res = 300)

library(circlize)
library(viridis)
circos.clear()
circos.par(start.degree = 90, track.margin = c(0.01, 0.01))
circos.initializeWithIdeogram(
  species = "hg38",
  ideogram.height = 0.1,        
  labels.cex = 1.0
)
my_colors <- magma(nrow(df), begin = 0, end = 0.8, alpha = 1) 
circos.genomicLabels(region1, 
                     labels = df[[GENE1]], 
                     side = "inside", 
                     cex = 0.5, 
                     col = my_colors,    
                     line_col = my_colors, 
                     connection_height = convert_height(5, "mm"))

circos.genomicLabels(region2, 
                     labels = df[[GENE2]], 
                     side = "outside", 
                     col = my_colors, 
                     line_col = my_colors, 
                     cex = 0.45)
circos.genomicLink(region1, region2, 
                   col = my_colors, 
                   border = NA)
circos.clear()
dev.off()
