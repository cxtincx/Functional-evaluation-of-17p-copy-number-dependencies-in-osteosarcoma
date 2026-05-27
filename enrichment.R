
library(AnnotationDbi)

library(DESeq2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db) 
library(dplyr)

go_genes <- df_HOS_OS1 %>%
  dplyr::filter(Group == "NRH-OS1")
go_final_df <- go_genes_interest %>%
  rownames_to_column(var = "Gene_Symbol")
my_gene_vector <- as.character(go_final_df$Gene_Symbol)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2) # For customizing the plot theme

gene_ids <- bitr(rownames(go_genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_region <- bitr(rownames(df_HOS_OS1), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_clean <- as.character(universe_region$ENTREZID)

go_results <- enrichGO(gene          = gene_ids$ENTREZID,
                       universe = universe_clean, 
                       OrgDb         = org.Hs.eg.db,
                       keyType       = "ENTREZID",
                       ont           = "BP", 
                       pAdjustMethod = "BH", 
                       pvalueCutoff  = 0.05)

go_results@result$Description <- str_to_sentence(go_results@result$Description)

library(clusterProfiler)
library(ggplot2)
library(stringr)

go_df <- as.data.frame(go_results)

parse_ratio <- function(ratio_str) {
  sapply(strsplit(ratio_str, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

go_df$FoldEnrichment <- parse_ratio(go_df$GeneRatio) / parse_ratio(go_df$BgRatio)

top_enriched <- head(go_df[go_df$p.adjust < 0.05, ], 10)
top_enriched <- top_enriched[order(top_enriched$FoldEnrichment, decreasing = FALSE), ]
top_enriched$Description <- str_to_sentence(top_enriched$Description)
top_enriched$Description <- factor(top_enriched$Description, levels = top_enriched$Description)

amp_genes_clean <- amp_genes$gene
amp_mapped <- bitr(amp_genes_clean, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)
amp_entrez <- amp_mapped$ENTREZID[!is.na(amp_mapped$ENTREZID)]

top_enriched$HasAmpGene <- sapply(top_enriched$geneID, function(x) {
  pathway_genes <- unlist(strsplit(x, "/"))
  any(pathway_genes %in% amp_entrez)
})

top_enriched$LabelText <- ifelse(top_enriched$HasAmpGene, "*", "")

p <- ggplot(top_enriched, aes(x = FoldEnrichment, y = Description, fill = p.adjust)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = LabelText), hjust = -0.2, vjust = 0.8, size = 10, color = "black") +
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") +
  labs(x = "Fold Enrichment") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 18, face = "plain", color = "black"),
    axis.text.x = element_text(size = 18, face = "plain"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16)
  )

print(p)
ggsave(filename = "/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/amp_gene_ontology_OS1_HOS.png", plot = p, width = 11, height = 5, dpi = 300)



amp_genes_clean <- as.character(amp_genes$gene)
amp_genes_clean <- trimws(amp_genes_clean)
amp_genes_clean <- amp_genes_clean[!is.na(amp_genes_clean) & amp_genes_clean != ""]

amp_mapped <- bitr(amp_genes_clean, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)

entrez_to_symbol <- amp_mapped$SYMBOL
names(entrez_to_symbol) <- amp_mapped$ENTREZID
amp_entrez <- amp_mapped$ENTREZID[!is.na(amp_mapped$ENTREZID)]

overlap_summary <- do.call(rbind, lapply(top_enriched$geneID, function(x) {
  pathway_genes <- unlist(strsplit(x, "/"))
  matching_ids <- intersect(pathway_genes, amp_entrez)
  matching_symbols <- entrez_to_symbol[matching_ids]
  
  data.frame(
    Count = length(matching_ids),
    Genes = paste(matching_symbols, collapse = ", "),
    stringsAsFactors = FALSE
  )
}))

pathway_gene_table <- data.frame(
  Pathway = top_enriched$Description,
  Amp_Gene_Count = overlap_summary$Count,
  Overlapping_Genes = overlap_summary$Genes,
  stringsAsFactors = FALSE
)

print(pathway_gene_table)

write.csv(pathway_gene_table, 
          file = "/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/amp_genes_pathway_summary.csv", 
          row.names = FALSE)



##DO 
# 1. Load Required Libraries
library(AnnotationDbi)
library(DESeq2)
library(clusterProfiler)
library(DOSE)               # Crucial for enrichDO
library(enrichplot)
library(org.Hs.eg.db) 
library(dplyr)
library(ggplot2)
library(stringr)
library(tibble)

go_genes <- df_HOS_OS1 %>%
  dplyr::filter(Group == "NRH-OS1")
gene_ids <- bitr(rownames(go_genes), 
                  fromType = "ENSEMBL", 
                  toType   = "ENTREZID", 
                  OrgDb    = org.Hs.eg.db)

do_results <- enrichDO(gene = gene_ids$ENTREZID,
                       ont = "HDO",       
                       pAdjustMethod = "BH", 
                       pvalueCutoff  = 0.05,
                       qvalueCutoff  = 0.2)        

amp_genes_clean <- amp_genes$gene
amp_mapped <- bitr(amp_genes_clean, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)

amp_entrez <- amp_mapped$ENTREZID[!is.na(amp_mapped$ENTREZID)]

top_do_enriched$HasAmpGene <- sapply(top_do_enriched$geneID, function(x) {
  pathway_genes <- unlist(strsplit(x, "/"))
  any(pathway_genes %in% amp_entrez)
  })
top_do_enriched$LabelText <- ifelse(top_do_enriched$HasAmpGene, "*", "")
p_do <- ggplot(top_do_enriched, aes(x = FoldEnrichment, y = Description, fill = p.adjust)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = LabelText), hjust = -0.2, vjust = 0.8, size = 10, color = "black") +
  scale_fill_gradient(low = "red", high = "blue", name = "p.adjust") +
  labs(x = "Fold Enrichment", title = "Disease Ontology Enrichment Analysis") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 14, face = "plain", color = "black"),
    axis.text.x = element_text(size = 18, face = "plain"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )

print(p_do)

ggsave(filename = "/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/disease_ontology.png", plot = p_do, width = 11, height = 5, dpi = 300)


amp_genes_clean <- as.character(amp_genes$gene)
amp_genes_clean <- trimws(amp_genes_clean)
amp_genes_clean <- amp_genes_clean[!is.na(amp_genes_clean) & amp_genes_clean != ""]

amp_mapped <- bitr(amp_genes_clean, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db, drop = FALSE)

entrez_to_symbol <- amp_mapped$SYMBOL
names(entrez_to_symbol) <- amp_mapped$ENTREZID
amp_entrez <- amp_mapped$ENTREZID[!is.na(amp_mapped$ENTREZID)]

overlap_summary <- do.call(rbind, lapply(top_do_enriched$geneID, function(x) {
  pathway_genes <- unlist(strsplit(x, "/"))
  matching_ids <- intersect(pathway_genes, amp_entrez)
  matching_symbols <- entrez_to_symbol[matching_ids]
  
  data.frame(
    Count = length(matching_ids),
    Genes = paste(matching_symbols, collapse = ", "),
    stringsAsFactors = FALSE
  )
}))

pathway_gene_table <- data.frame(
  Pathway = top_do_enriched$Description,
  Amp_Gene_Count = overlap_summary$Count,
  Overlapping_Genes = overlap_summary$Genes,
  stringsAsFactors = FALSE
)

print(pathway_gene_table)

