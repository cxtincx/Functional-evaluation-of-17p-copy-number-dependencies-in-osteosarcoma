
library("tximport")
library("dplyr")

## assign your working directory
setwd ("/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1")
## List all directories containing data
samples <- list.files(path = "./Salmon", full.names = T, pattern = "Salmon_.*\\.sf$")
## Obtain a vector of all filenames including the path
files <- samples
## list all the files in the console
files
##assign a shorter name for each element
names(files) <- c(
  "HOS_1",
  "HOS_2",
  "HOS_3",
  "OS1_1",
  "OS1_2",
  "OS1_3"
)
## Read the annotation file
# Install if needed
# install.packages("BiocManager")

library(rtracklayer)

# Path to the GENCODE GTF you downloaded
gtf_file <- "/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/gencode.v46.basic.annotation.gtf"

# Import GTF
gtf <- rtracklayer::import(gtf_file)

# Keep transcript entries
tx <- gtf[gtf$type == "transcript"]

# Extract IDs
tx_id   <- mcols(tx)$transcript_id
gene_id <- mcols(tx)$gene_id

# Strip version numbers (e.g. ENST00000.1 -> ENST00000)
strip_version <- function(x) sub("\\.[0-9]+$", "", x)

tx2gene <- unique(
  data.frame(
    TXNAME = strip_version(as.character(tx_id)),
    GENEID = strip_version(as.character(gene_id)),
    stringsAsFactors = FALSE
  )
)

# Save for tximport
write.csv(tx2gene, "tx2gene_gencode_grch38.csv",
          row.names = FALSE, quote = FALSE)

## if you want to view this file type head (tx2gene)
# Run tximport
txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("TXNAME", "GENEID")], countsFromAbundance="lengthScaledTPM",
                ignoreTxVersion = TRUE)
## Check the output
head(txi[["counts"]])
## Extract the counts, round the values and change the list to a data frame
counts <- txi$counts %>% round() %>% data.frame()
## Change the column names
##colnames (counts) <- c("5aza_rep1", "5aza_rep2", "DMSO_rep1", "DMSO_rep2" )
##View(counts)
### Save the results for the next session using the write.table command
write.table(counts, "/Users/catincaapostol/Desktop/Master's/rnaseq nrhos1/Salmon_counts_tximport.txt", sep ="\t", quote = F)
####keep a note of the packages and their versions
##sessionInfo()
##DGE CODE 

#Read the count data generated in the previous code
#Load the Sample information (metadata).
counts <- read.table(file = "/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/Salmon_counts_tximport.txt", header= TRUE, check.names = FALSE)
metadata <- data.frame(
  condition = c("HOS", "HOS", "HOS",
                "NRH-OS1", "NRH-OS1", "NRH-OS1"),
  row.names = colnames(counts)
)

metadata$condition <- factor(
  metadata$condition,
  levels = c("HOS", "NRH-OS1")
)
metadata$replicate <- c(1, 2, 3, 1, 2, 3)

##check if the sample name is matching the metadata
all(rownames(metadata)%in% colnames(counts))
all(rownames(metadata)== colnames(counts))
#If the order of rows and columns is not the same, try do the following
##counts<- counts[, row.names(metadata)]

design <- as.formula(~condition)
model<- model.matrix(design, data= metadata)
keep <- rowSums(counts)>5
countdata<- counts[keep,]
countdata<- as.matrix(countdata)
#### Create DESeq2Dataset object
library(DESeq2)
dds.raw <- DESeqDataSetFromMatrix(countData = countdata,
                                  colData = metadata,
                                  design = design)
##Perform the differential gene expression
library(apeglm)
dds <- DESeq(dds.raw)
res_05<- results(dds,alpha= 0.05)
resLFC_05 <- lfcShrink(dds, coef="condition_NRH.OS1_vs_HOS", type="apeglm", res= res_05)

DEG <- as.data.frame(resLFC_05) # assign the results to data frame DEG
# assign symbol(common names) and ENTREZ id according to ENSEMBL id
library(AnnotationDbi)
library(org.Hs.eg.db)
DEG$symbol<- mapIds(org.Hs.eg.db,
                    keys= rownames(DEG),
                    column = "SYMBOL",
                    keytype = "ENSEMBL",
                    multiVals = "first")
DEG$entrez<- mapIds(org.Hs.eg.db,
                    keys= rownames(DEG),
                    column = "ENTREZID",
                    keytype = "ENSEMBL",
                    mutiVals= "first")
# remove genes that don't have a common name and those with duplicated gene name
DEG_symbol<- DEG[is.na(DEG$symbol)== FALSE,]
dim(DEG_symbol)
DEG_symbol<- DEG_symbol[!duplicated(DEG_symbol$symbol),]

#volcano plot 
library(ggplot2)
library(ggrepel)   # optional, for nice labels

# Use the shrunken results if you have them
# resLFC_05 should be a DESeq2 results-like object or data.frame
df_HOS_OS1 <- as.data.frame(DEG_symbol)
names(df_HOS_OS1)[names(df_HOS_OS1) == "group"] <- "Group"
# make sure it has the columns: log2FoldChange and padj (or pvalue)
head(df_HOS_OS1)

# Add plotting columns (handle NAs)
df_HOS_OS1 <- df_HOS_OS1[!is.na(df_HOS_OS1$padj), ]
df_HOS_OS1$negLog10Padj <- -log10(df_HOS_OS1$padj)

# Decide thresholds
padj_cutoff <- 0.05

# classification for coloring
df_HOS_OS1$Group <- "NS"  # not significant
df_HOS_OS1$Group[ df_HOS_OS1$padj < padj_cutoff & df_HOS_OS1$log2FoldChange > 0 ] <- "NRH-OS1"
df_HOS_OS1$Group[ df_HOS_OS1$padj < padj_cutoff & df_HOS_OS1$log2FoldChange < 0 ] <- "HOS"

# keep gene names for label examples (row names)
df_HOS_OS1$genes <- rownames(df_HOS_OS1)

# Basic volcano
p <- ggplot(df_HOS_OS1, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = Group), alpha = 1, size = 3) +
  scale_color_manual(values = c("HOS" = "steelblue", "NRH-OS1" = "firebrick", "NS" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey60") +
  labs(
    x = "log2 Fold Change (HOS / NRH-OS1)",
    y = "-log10(adj. p-value)",
    title = "Differential gene expression: HOS vs NRH-OS1",
    subtitle = paste0("Left = HOS ", " Right = NRH-OS1")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16), 
    axis.text.x = element_text(size = 16), 
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 16)
  )

top_genes <- dplyr::filter(df_HOS_OS1, Group != "NS") %>%
  dplyr::arrange(desc(negLog10Padj)) %>%
  head(40)
genes_to_add <- c("MYC", "TP53", "UBA6", "NUB1", "UBB", "UBC")
top_genes <- dplyr::bind_rows(
  top_genes,
  dplyr::filter(df_HOS_OS1, symbol %in% genes_to_add)
) %>%
  dplyr::distinct(symbol, .keep_all = TRUE)

p <- p + geom_text_repel(
  data = top_genes,
  aes(label = symbol),
  size = 4.5,
  max.overlaps = 5
)

print(p)
ggsave("/Users/catincaapostol/Desktop/Master's/project/volcano_OS1_OHS.png", p, width = 7, height = 4.5, dpi = 300)


#chr17 zoom in ------------------------------------------------------------------------------------------
###### HOS VS OS1 chr17p12-11.2


# 1) Clean annotation_df: remove NA/blank gene names and exact duplicate rows
library(dplyr)

# 0) Quick diagnostics (run this if you want to inspect names)
print(names(ann_df))

# 1) Normalize column names (trim whitespace) — helps if column names include stray spaces
names(ann_df) <- trimws(names(ann_df))

# 2) Clean annotation_df: remove NA/blank gene names and exact duplicate rows
annotation_clean_HOS_OS1 <- ann_df %>%
  filter(!is.na(`Gene name`) & trimws(`Gene name`) != "") %>%
  distinct()

# 3) Make sure the key column exists, then collapse other columns (if any)
if (!("Gene name" %in% names(annotation_clean_HOS_OS1))) {
  stop("Column 'Gene name' not found in annotation_df. Run names(annotation_df) to inspect column names.")
}

other_cols <- setdiff(names(annotation_clean_HOS_OS1), "Gene name")

if (length(other_cols) > 0) {
  annotation_collapsed_HOS_OS1 <- annotation_clean_HOS_OS1 %>%
    group_by(`Gene name`) %>%
    summarise(
      across(all_of(other_cols),
             ~ {
               vals <- unique(na.omit(as.character(.)))
               if (length(vals) == 0) "" else paste(vals, collapse = "; ")
             },
             .names = "{.col}"),
      .groups = "drop"
    )
} else {
  # only Gene name column present
  annotation_collapsed_HOS_OS1 <- annotation_clean_HOS_OS1 %>% distinct()
}


#VOLCANO PLOT 
library(ggplot2)
library(ggrepel)   # optional, for nice labels

# Use the shrunken results if you have them
# resLFC_05 should be a DESeq2 results-like object or data.frame
df_location_HOS_OS1 <- df_HOS_OS1
df_location_HOS_OS1 <- as.data.frame(DEG_symbol)

# make sure it has the columns: log2FoldChange and padj (or pvalue)
head(df_location_HOS_OS1)

# Add plotting columns (handle NAs)
df_location_HOS_OS1 <- df_location_HOS_OS1[!is.na(df_location_HOS_OS1$padj), ]
df_location_HOS_OS1$negLog10Padj <- -log10(df_location_HOS_OS1$padj)

# Decide thresholds
padj_cutoff <- 0.05

# classification for coloring
df_location_HOS_OS1$group <- "NS"  # not significant
df_location_HOS_OS1$group[ df_location_HOS_OS1$padj < padj_cutoff & df_location_HOS_OS1$log2FoldChange > 0 ] <- "OHS"
df_location_HOS_OS1$group[ df_location_HOS_OS1$padj < padj_cutoff & df_location_HOS_OS1$log2FoldChange < 0 ] <- "HOS"

# keep gene names for label examples (row names)
df_location_HOS_OS1$genes <- rownames(df_location_HOS_OS1)

# 4) Merge into df (left join to keep only df rows but attach annotation where available)
df_location_HOS_OS1 <- df_location_HOS_OS1 %>%
  left_join(annotation_collapsed, by = c("symbol" = "Gene name")) %>%
  distinct()  # final deduplication of entire rows
# Basic volcano
p <- ggplot(df_location_HOS_OS1, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = group), alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("HOS" = "steelblue", "OS1" = "firebrick", "NS" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey60") +
  labs(
    x = "log2 Fold Change (HOS / OS1)",
    y = "-log10(adj. p-value)",
    title = "Volcano plot: HOS vs OS1",
    subtitle = paste0("Left = HOS Right = OS1")
  ) +
  theme_minimal(base_size = 14)

# Optionally label top genes (highest -log10(padj) or highest abs(LFC))
top_genes_HOS_OS1 <- df_location_HOS_OS1 %>%
  dplyr::filter(
    group != "NS",
    `Chromosome/scaffold name` == 17,
    `Karyotype band` %in% c("p11.2", "p12")
  ) %>%
  dplyr::arrange(desc(negLog10Padj))
  

p <- p + geom_text_repel(
  data = top_genes_HOS_OHS,
  aes(label = symbol),
  size = 3,
  max.overlaps = 25
)
print(p)
ggsave("volcano_HOS_OHS_17p.png", p, width = 16, height = 6, dpi = 300)

#keep OHS and OS1 upreg 
library(dplyr)

# symbols up in OHS (in the HOS_OHS comparison) and symbols HOS in the other df
symbols_up_in_OHS <- df_location_HOS_OHS %>%
  dplyr::filter(group == "OHS") %>%
  dplyr::pull(symbol) %>%
  unique()

symbols_HOS_in_both <- df_location_HOS_OHS %>%
  dplyr::filter(group == "HOS") %>%
  dplyr::pull(symbol) %>%
  unique()

# keep only rows from df_location_HOS (all its columns) that meet either condition:
#  - group == "OS1" in df_location_HOS AND symbol is up in OHS in df_location_HOS_OHS
#  - group == "HOS" in df_location_HOS AND symbol is HOS in df_location_HOS_OHS
genes_keep_df <- df_location_HOS %>%
  dplyr::filter(
    (group == "OS1" & symbol %in% symbols_up_in_OHS) |
      (group == "HOS" & symbol %in% symbols_HOS_in_both)
  )

normalised_counts <- txi$counts
