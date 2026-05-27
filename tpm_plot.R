library("tximport")
library("dplyr")
library(scales)
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

txi <- tximport(files, type="salmon", tx2gene=tx2gene[,c("TXNAME", "GENEID")], countsFromAbundance="lengthScaledTPM",
                ignoreTxVersion = TRUE)

head(txi[["counts"]])
counts <- txi$counts %>% round() %>% data.frame()

#meank the tpms 
tpm_hos_df <- mutate(counts, HOS_TPM = rowMeans(dplyr::select(counts, c("HOS_1", "HOS_2", "HOS_3")), na.rm = TRUE)) %>% round() %>% data.frame()

tpm_df <- mutate(tpm_hos_df, OS1_TPM = rowMeans(dplyr::select(counts, c("OS1_1", "OS1_2", "OS1_3")), na.rm = TRUE)) %>% round() %>% data.frame()

#replace gene name 

library(AnnotationDbi)
library(org.Hs.eg.db)
tpm_genes <- tpm_df
tpm_genes$gene_name <- mapIds(
  org.Hs.eg.db,
  keys      = rownames(tpm_df),   # Accesses the first column regardless of its name
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# remove genes that don't have a common name and those with duplicated gene name
tpm_final<- tpm_genes[is.na(tpm_genes$gene_name)== FALSE,]
dim(tpm_final)
tpm_final<- tpm_final[!duplicated(tpm_final$gene_name),]

#make one df 
wgs_tpm_df_common <- tpm_final %>%
  dplyr::select(
    HOS_TPM,
    OS1_TPM,
    gene_name
  ) %>%
  dplyr::inner_join(
    OS1_copy_number_outrem %>%
      dplyr::select(gene, cn, order),
    by = c("gene_name" = "gene")
  )

p <- ggplot(
  wgs_tpm_df_common,
  aes(x = cn, y = OS1_TPM)
) +
  geom_point(size = 3, alpha = 1, color = "grey") +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_x_continuous(
    name = "Copy number (CN)",
    breaks = seq(0, 25, by = 5)
  ) +
  scale_y_continuous(trans = pseudo_log_trans(base = 10, sigma = 25), name = "TPM"
  ) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.line.y = element_line(color = "black", size = 0.8),
    axis.line.x = element_line(color = "black", size = 0.8),
    axis.ticks = element_line(color = "black"),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20)
  )

print(p)

library(ggpmisc)
library(scales)
b <- 10
s <- 1
formula <- asinh(y / (2 * s)) / log(b) ~ x

p_r <- p +
  stat_poly_line(
    formula = formula, 
    color = "black",
    linetype = "dashed",
    linewidth = 0.5, 
    se = FALSE            
  ) +
  stat_poly_eq(
    formula = formula,
    aes(label = after_stat(rr.label)), 
    label.x = "right", 
    label.y = "top", 
    parse = TRUE,
    size = 4.5
  )

print(p_r)

#filter for 17p
wgs_tpm_df_common_17 <- tpm_final %>%
  dplyr::select(
    HOS_TPM,
    OS1_TPM,
    gene_name
  ) %>%
  dplyr::inner_join(
    chr17_OS1_copy_number_outrem %>%
      dplyr::select(gene, cn, order_17, location),
    by = c("gene_name" = "gene")
  )
wgs_tpm_df_common_17 <- wgs_tpm_df_common_17[wgs_tpm_df_common_17$gene_name %in% amp_genes$gene,]

#add points


p_final <- p_r + 
  geom_point(data = wgs_tpm_df_common_17, aes(x = cn, y = OS1_TPM), size = 3, alpha = 1, color = "firebrick") +
  geom_text_repel(
    data = wgs_tpm_df_common_17,      
    aes(x = cn, y = OS1_TPM, label = gene_name),
    size = 5,
    max.overlaps = 6,
    show.legend = FALSE, 
    color = "black"
  ) 


print(p_final)

p_final_r <- p_final +
  stat_poly_line(
    data = wgs_tpm_df_common_17, 
    formula = formula, 
    color = "firebrick",   
    linetype = "dashed",  
    linewidth = 0.5, 
    se = FALSE            
  ) +
  stat_poly_eq(
    data = wgs_tpm_df_common_17,
    formula = formula,
    aes(label = after_stat(rr.label)), 
    label.x = "right", 
    label.y = 0.85,
    color = "firebrick",
    parse = TRUE,
    size = 4.5
  )

print(p_final_r)
ggsave("/Users/catincaapostol/Desktop/Master's/project/os1 analysis plots /17p_OS1_TPM_CN.png", p_final_r, width = 7, height = 5, dpi = 300)

