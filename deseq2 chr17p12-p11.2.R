library("tximport")
library("dplyr")

setwd ("/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1")
samples <- list.files(path = "./Salmon", full.names = T, pattern = "Salmon_.*\\.sf$")
files <- samples
files
names(files) <- c(
  "HOS_1",
  "HOS_2",
  "HOS_3",
  "OS1_1",
  "OS1_2",
  "OS1_3"
)

library(rtracklayer)
gtf_file <- "/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/gencode.v46.basic.annotation.gtf"


gtf <- rtracklayer::import(gtf_file)


tx <- gtf[gtf$type == "transcript"]

tx_id   <- mcols(tx)$transcript_id
gene_id <- mcols(tx)$gene_id

strip_version <- function(x) sub("\\.[0-9]+$", "", x)

tx2gene <- unique(
  data.frame(
    TXNAME = strip_version(as.character(tx_id)),
    GENEID = strip_version(as.character(gene_id)),
    stringsAsFactors = FALSE
  )
)

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
write.table(counts, "/Users/catincaapostol/Desktop/Master's/project/rnaseq nrhos1/Salmon_counts_tximport.txt", sep ="\t", quote = F)
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

DEG <- as.data.frame(resLFC_05) 
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
DEG_symbol<- DEG[is.na(DEG$symbol)== FALSE,]
dim(DEG_symbol)
DEG_symbol<- DEG_symbol[!duplicated(DEG_symbol$symbol),]


padj_cutoff <- 0.05
df_HOS_OS1 <- as.data.frame(DEG_symbol)
df_HOS_OS1$Group <- "NS"
df_HOS_OS1$Group[ df_HOS_OS1$padj < padj_cutoff & df_HOS_OS1$log2FoldChange > 0 ] <- "NRH-OS1"
df_HOS_OS1$Group[ df_HOS_OS1$padj < padj_cutoff & df_HOS_OS1$log2FoldChange < 0 ] <- "HOS"

df_HOS_OS1 <- df_HOS_OS1[!is.na(df_HOS_OS1$padj), ]
df_HOS_OS1$negLog10Padj <- -log10(df_HOS_OS1$padj)
head(df_HOS_OS1)

df_HOS_OS1_amp <- df_HOS_OS1 %>%
  dplyr::filter(symbol %in% amp_genes$gene)

#volcano plot 
library(ggplot2)
library(ggrepel)
df_HOS_OS1_17p <- df_HOS_OS1_17p[!is.na(df_HOS_OS1_17p$padj), ]
df_HOS_OS1_17p$negLog10Padj <- -log10(df_HOS_OS1_17p$padj)

head(df_HOS_OS1_17p)


# keep gene names for label examples (row names)
df_HOS_OS1_amp$symbol <- rownames(df_HOS_OS1_amp)
head(df_HOS_OS1_amp)
# Basic volcano
p <- ggplot(df_HOS_OS1_amp, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = Group), alpha = 1, size = 3) +
  scale_color_manual(values = c("HOS" = "steelblue", "NRH-OS1" = "firebrick", "NS" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey60") +
  labs(
    x = "log2 Fold Change (HOS / OS1)",
    y = "-log10(adj. p-value)",
    title = "Volcano plot: HOS vs OS1",
    subtitle = paste0("Left = HOS ", " Right = NRH-OS1")
  ) +
  theme_minimal(base_size = 14) + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16), 
    axis.text.x = element_text(size = 16), 
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 16),
    plot.subtitle = element_text(size = 20)
  )


# Optionally label top genes (highest -log10(padj) or highest abs(LFC))
p <- p + geom_text_repel(
  data = df_HOS_OS1_amp,
  aes(label = symbol),
  size = 4.5,
  max.overlaps = 6
)

#label all genes 


# print and/or save
print(p)
ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 analysis plots /volcano_OS1_HOS_amp.png", p, width = 7, height = 4.5, dpi = 300)


##DGE after filtering 
keep_genes <- rownames(countdata) %in% amp_genes$ensembl_gene_id
countdata_amp_genes <- countdata[keep_genes, ]

dds.raw <- DESeqDataSetFromMatrix(countData = countdata_amp_genes,
                                  colData = metadata,
                                  design = design)
##Perform the differential gene expression
library(apeglm)
dds <- DESeq(dds.raw)
res_05<- results(dds,alpha= 0.05)
resLFC_05 <- lfcShrink(dds, coef="condition_NRH.OS1_vs_HOS", type="apeglm", res= res_05)

DEG <- as.data.frame(resLFC_05) 
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
DEG_symbol<- DEG[is.na(DEG$symbol)== FALSE,]
dim(DEG_symbol)
DEG_symbol<- DEG_symbol[!duplicated(DEG_symbol$symbol),]


padj_cutoff <- 0.05
df_HOS_OS1 <- as.data.frame(DEG_symbol)
df_HOS_OS1$Group <- "NS"
df_HOS_OS1$Group[ df_HOS_OS1$padj < padj_cutoff & df_HOS_OS1$log2FoldChange > 0 ] <- "NRH-OS1"
df_HOS_OS1$Group[ df_HOS_OS1$padj < padj_cutoff & df_HOS_OS1$log2FoldChange < 0 ] <- "HOS"

df_HOS_OS1 <- df_HOS_OS1[!is.na(df_HOS_OS1$padj), ]
df_HOS_OS1$negLog10Padj <- -log10(df_HOS_OS1$padj)
head(df_HOS_OS1)

df_HOS_OS1_amp <- df_HOS_OS1 %>%
  dplyr::filter(symbol %in% amp_genes$gene)

#volcano plot 
library(ggplot2)
library(ggrepel)
df_HOS_OS1_17p <- df_HOS_OS1_17p[!is.na(df_HOS_OS1_17p$padj), ]
df_HOS_OS1_17p$negLog10Padj <- -log10(df_HOS_OS1_17p$padj)

head(df_HOS_OS1_17p)


# keep gene names for label examples (row names)
df_HOS_OS1_amp$symbol <- rownames(df_HOS_OS1_amp)
head(df_HOS_OS1_amp)
# Basic volcano
p <- ggplot(df_HOS_OS1_amp, aes(x = log2FoldChange, y = negLog10Padj)) +
  geom_point(aes(color = Group), alpha = 1, size = 3) +
  scale_color_manual(values = c("HOS" = "steelblue", "NRH-OS1" = "firebrick", "NS" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey60") +
  labs(
    x = "log2Fold-Change",
    y = "-log10(adj. p-value)",
    title = "Volcano plot: HOS vs OS1",
    subtitle = paste0("Left = HOS ", " Right = NRH-OS1")
  ) +
  theme_minimal(base_size = 14) + 
  scale_x_continuous(limits = c(-10, 10), breaks = c(-10, -7.5, -5, -2.5, 2.5, 5, 7.5, 10)) + 
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y = element_text(size = 16), 
    axis.text.x = element_text(size = 16), 
    legend.title = element_text(size = 20), 
    legend.text = element_text(size = 16),
    plot.subtitle = element_text(size = 20)
  )


# Optionally label top genes (highest -log10(padj) or highest abs(LFC))
p <- p + geom_text_repel(
  data = df_HOS_OS1_amp,
  aes(label = symbol),
  size = 4.5,
  max.overlaps = 7
)

#label all genes 


# print and/or save
print(p)

ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 analysis plots /volcano_OS1_HOS_amp.png", p, width = 7, height = 5, dpi = 300)

masters_theme <- p$theme
