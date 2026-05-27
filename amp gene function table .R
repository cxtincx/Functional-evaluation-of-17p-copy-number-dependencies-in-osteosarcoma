library(org.Hs.eg.db)
library(AnnotationDbi)

amp_genes_function <- select(org.Hs.eg.db, 
       keys = amp_genes$gene, 
       columns = c("GENENAME", "GO"), 
       keytype = "SYMBOL")
amp_genes_function <- amp_genes_function[!duplicated(amp_genes_function$SYMBOL), ]

library(dplyr)
library(stringr)

amp_genes_function <- amp_genes_function %>%
  mutate(GENENAME = str_replace(GENENAME, "^([a-zA-Z])", toupper))

names(amp_genes_function)[names(amp_genes_function) == "GENENAME"] <- "Gene Role"

amp_genes_function <- amp_genes_function %>%
  dplyr::select("SYMBOL", "Gene Role")
rownames(amp_genes_function) <- NULL

write.csv(amp_genes_function, "/Users/catincaapostol/Desktop/Master's/project/wgs os1/Amplified genes function.csv")

