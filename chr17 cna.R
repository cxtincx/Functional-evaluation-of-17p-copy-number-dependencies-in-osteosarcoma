#list of protein coding genes
Protein_type<- read.csv("/Users/catincaapostol/Desktop/Master's/project/mart_export_old.txt")
names(Protein_type)[names(Protein_type) == "Gene.stable.ID"] <- "gene_id"
names(Protein_type)[names(Protein_type) == "Gene.name"] <- "gene_name"
Protein_type_coding <- subset(Protein_type, Gene.type == "protein_coding")
Protein_type_coding [Protein_type_coding  == ""] <- NA
Protein_type_coding <- Protein_type_coding[!is.na(Protein_type_coding$gene_name), ]
protein_coding_list<-Protein_type_coding$gene_name
names(Protein_type_coding)[names(Protein_type_coding) == "gene_name"] <- "gene"
Protein_type_coding<-Protein_type_coding[Protein_type_coding$Chromosome.scaffold.name %in% c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22"),]
Protein_type_coding<-Protein_type_coding[,colnames(Protein_type_coding) %in% c("gene","Gene.start..bp.", "Gene.end..bp.","Chromosome.scaffold.name")]

#Filter to include only protein coding genes
OS1_CNV_calling_coding<-OS1_CNV_calling

OS1_CNV_calling_coding$coding <- strsplit(OS1_CNV_calling_coding$gene, ",")
for (i in seq_len(nrow(OS1_CNV_calling_coding))) {
  # Filter genes in the current row that are in protein_coding_list
  OS1_CNV_calling_coding$coding[[i]] <- OS1_CNV_calling_coding$coding[[i]][OS1_CNV_calling_coding$coding[[i]] %in% protein_coding_list]
}
OS1_CNV_calling_coding_filtered <- OS1_CNV_calling_coding[sapply(OS1_CNV_calling_coding$coding, length) > 0, ]

OS1_CNV_calling_coding_filtered$coding_new <- sapply(OS1_CNV_calling_coding_filtered$coding, function(x) paste(x, collapse = ","))


#cell in the coding_new column contains multiple values separated by commas (,), each value will be split into a new row. 
OS1_CNV_calling_coding_filtered_unique <- OS1_CNV_calling_coding_filtered %>%
  separate_rows(coding_new, sep = ",")
#Deleted duplicated genes with the same cn 
OS1_CNV_calling_coding_filtered_unique_filtered <- OS1_CNV_calling_coding_filtered_unique %>%
  distinct(coding_new, cn, .keep_all = TRUE)

#Some genes have two different copy number values - this puts them in a table
OS1_repeated <- OS1_CNV_calling_coding_filtered_unique_filtered %>%
  group_by(coding_new) %>%
  filter(n() > 1) %>%
  ungroup()


#create a list of the repeated genes (genes with more than one value for copy number)
OS1_repeated_genes<- OS1_repeated$coding_new
OS1_repeated_genes<- unique(OS1_repeated_genes)

#Determine the mean cn for the genes with more than one cn value, if it is a decimal CN, round up
cn_OS1_mean <- aggregate(cn ~ coding_new, data = OS1_repeated, FUN = mean)
cn_OS1_mean$cn_round <- ceiling(cn_OS1_mean$cn)

#Create a new df to extract the information for each gene from OHS_repeated and format the df to match OHS_coding_copynumber_fil
OS1_repeated_unique <- OS1_repeated[!duplicated(OS1_repeated$coding_new), ]
OS1_repeated_unique <- OS1_repeated_unique[, !colnames(OS1_repeated_unique) %in% c("gene", "coding", "cn")]
OS1_meancn<- merge(OS1_repeated_unique, cn_OS1_mean, by = "coding_new")

OS1_meancn <- OS1_meancn[, !colnames(OS1_meancn) %in% c("cn")]
names(OS1_meancn)[names(OS1_meancn) == "coding_new"] <- "gene"
names(OS1_meancn)[names(OS1_meancn) == "cn_round"] <- "cn"
OS1_meancn<-OS1_meancn[,c("chromosome","start", "end", "log2", "ci_hi", "ci_lo", "cn", "depth", "probes", "weight", "gene")]


#delete rows of genes with 2 different cn values form the original df 
OS1_CNV_calling_removedup<-OS1_CNV_calling_coding_filtered_unique_filtered
OS1_CNV_calling_removedup<-OS1_CNV_calling_removedup[!OS1_CNV_calling_removedup$coding_new %in% OS1_repeated_genes,]
OS1_CNV_calling_removedup$coding<-NULL
OS1_CNV_calling_removedup$gene<-NULL
names(OS1_CNV_calling_removedup)[names(OS1_CNV_calling_removedup) == "coding_new"] <- "gene"

OS1_final_cn<-rbind(OS1_CNV_calling_removedup,OS1_meancn )

#Ordering the final df

OS1_final_cn$chromosome <- gsub("chr", "", OS1_final_cn$chromosome)
OS1_final_cn<-OS1_final_cn[!OS1_final_cn$chromosome =="X",]
OS1_final_cn<-OS1_final_cn[!OS1_final_cn$chromosome =="Y",]
OS1_final_cn$chromosome<-as.numeric(OS1_final_cn$chromosome)

OS1_final_cn_order<-merge(OS1_final_cn, Protein_type_coding, by = "gene")


setdiff(OS1_final_cn$gene,OS1_final_cn_order$gene)
OS1_final_cn_order <- OS1_final_cn_order[order(OS1_final_cn_order$Gene.start..bp.), ]
OS1_final_cn_order <- OS1_final_cn_order[order(OS1_final_cn_order$chromosome), ]

OS1_final_cn_order$order <- seq_along(OS1_final_cn_order$gene)

write.csv(OS1_final_cn_order,"/Users/catincaapostol/Desktop/Master's/project/wgs os1/OS1_cnploidy.csv")

OS1_copy_number<-OS1_final_cn_order[,colnames(OS1_final_cn_order) %in% c("gene", "cn","chromosome", "order")]

write.csv(OS1_copy_number,"/Users/catincaapostol/Desktop/Master's/project/wgs os1/OS1_processed_cn.csv")



smallest_order_OS1 <- tapply(OS1_copy_number$order, OS1_copy_number$chromosome, min)
largest_order_OS1 <- tapply(OS1_copy_number$order, OS1_copy_number$chromosome, max)

positions_OS1 <- data.frame(
  chromosome = names(smallest_order_OS1),
  smallest_order = smallest_order_OS1,
  largest_order = largest_order_OS1
)
positions_OS1$smallest_order<- positions_OS1$smallest_order -1 
positions_OS1$mid<- (positions_OS1$smallest_order + positions_OS1$largest_order) / 2

######
##Zoom in on chr17
Mart<-read.csv("/Users/catincaapostol/Desktop/Master's/project/hgnc_complete_set_subset_order_CHROMOSOME_kb.csv")
Mart_chr17<-Mart[Mart$Chromosome.scaffold.name == "17",]
Mart_chr17<-unique(Mart_chr17)
Mart_chr17<-Mart_chr17[!Mart_chr17$Gene.name =="",]
colnames(Mart_chr17)[colnames(Mart_chr17) == "Gene.name"] <- "gene"
cyto <- read.table(
  "/Users/catincaapostol/Desktop/Master's/project/wgs os1/cytoBand.txt",
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  fill = TRUE,
  quote = ""
)

colnames(cyto) <- c("chrom", "chromStart", "chromEnd", "band", "stain")
cyto <- cyto[grepl("^chr[0-9XY]+$", cyto$chrom), ]
cyto17 <- cyto[cyto$chrom == "chr17", ]


##OS1
chr17_OS1_copy_number_outrem<-OS1_copy_number_outrem[OS1_copy_number_outrem$chromosome == "17",]
chr17_OS1_copy_number_outrem<-merge(chr17_OS1_copy_number_outrem, Mart_chr17, by="gene")
## OS1
chr17_OS1_copy_number_outrem <- OS1_copy_number_outrem[
  OS1_copy_number_outrem$chromosome == "17", ]

chr17_OS1_copy_number_outrem <- merge(
  chr17_OS1_copy_number_outrem,
  Mart_chr17,
  by = "gene"
)

# order genes along chr17
chr17_OS1_copy_number_outrem <- chr17_OS1_copy_number_outrem[
  order(chr17_OS1_copy_number_outrem$order), ]

chr17_OS1_copy_number_outrem$order_17 <- seq_len(nrow(chr17_OS1_copy_number_outrem))
chr17_OS1_copy_number_outrem$order_17copy <- chr17_OS1_copy_number_outrem$order_17

possible_starts <- c(
  "Gene.start..bp.", "Gene.start..bp", "Gene.start", "start", "Start",
  "start.x", "start.y", "gene_start", "gene_start_bp", "location.start"
)
possible_ends <- c(
  "Gene.end..bp.", "Gene.end..bp", "Gene.end", "end", "End",
  "end.x", "end.y", "gene_end", "gene_end_bp", "location.end"
)

cols <- colnames(chr17_OS1_copy_number_outrem)

start_col <- intersect(possible_starts, cols)
end_col   <- intersect(possible_ends, cols)

start_col <- if(length(start_col) > 0) start_col[1] else NA_character_
end_col   <- if(length(end_col)   > 0) end_col[1]   else NA_character_

if (is.na(start_col)) {
  stop(
    "No gene start column found. Available columns:\n  ",
    paste(cols, collapse = ", "),
    "\n\nAdd a start column (e.g. Gene.start..bp.) or rename your coordinate columns."
  )
}

gene_start_num <- as.numeric(chr17_OS1_copy_number_outrem[[start_col]])
if (!is.na(end_col)) {
  gene_end_num <- as.numeric(chr17_OS1_copy_number_outrem[[end_col]])
  gene_mid <- ifelse(is.na(gene_end_num), gene_start_num, (gene_start_num + gene_end_num) / 2)
} else {
  gene_mid <- gene_start_num
}

if (all(is.na(gene_mid))) stop("Gene coordinate column found but contains only NA/NA-like values.")
if (length(gene_mid) != nrow(chr17_OS1_copy_number_outrem)) stop("Length mismatch between gene_mid and rows of chr17_OS1_copy_number_outrem.")

if (!exists("cyto17")) stop("cyto17 not found. Create cyto17 <- cyto[cyto$chrom == 'chr17', ] before this block.")
cyto17$chromStart <- as.numeric(cyto17$chromStart)
cyto17$chromEnd   <- as.numeric(cyto17$chromEnd)
chr17_OS1_copy_number_outrem$Karyotype.band <- vapply(
  gene_mid,
  FUN = function(pos) {
    if (is.na(pos)) return(NA_character_)
    idx <- which(pos >= cyto17$chromStart & pos <= cyto17$chromEnd)
    if (length(idx) == 0) return(NA_character_) else return(as.character(cyto17$band[idx[1]]))
  },
  FUN.VALUE = character(1)
)
cat("Karyotype.band assignment: rows =", nrow(chr17_OS1_copy_number_outrem), 
    "-> assigned bands =", sum(!is.na(chr17_OS1_copy_number_outrem$Karyotype.band)), 
    " / ", nrow(chr17_OS1_copy_number_outrem), "\n")
cat("Top bands (counts):\n")
print(head(sort(table(chr17_OS1_copy_number_outrem$Karyotype.band), decreasing = TRUE), 20))




# assign each gene to a cytoband based on order
cyto17$order_band_start <- seq_len(nrow(cyto17))


gene_mid <- (chr17_OS1_copy_number_outrem$start +
               chr17_OS1_copy_number_outrem$end) / 2

# ensure cytobands are ordered
cyto17 <- cyto17[order(cyto17$chromStart), ]
gene_mid <- (chr17_OS1_copy_number_outrem$Gene.start..bp. +
               chr17_OS1_copy_number_outrem$Gene.end..bp.) / 2

# map each gene midpoint to a cytoband
chr17_OS1_copy_number_outrem$Karyotype.band <- sapply(
  gene_mid,
  function(pos) {
    idx <- which(pos >= cyto17$chromStart & pos <= cyto17$chromEnd)
    if (length(idx) == 0) NA else cyto17$band[idx[1]]
  }
)



chr17_OS1_copy_number_outrem$order_17 <-
  ifelse(chr17_OS1_copy_number_outrem$order_17copy > 323,
         chr17_OS1_copy_number_outrem$order_17copy + 50,
         chr17_OS1_copy_number_outrem$order_17copy)

smallest_order_kb_OS1 <-
  tapply(chr17_OS1_copy_number_outrem$order_17,
         chr17_OS1_copy_number_outrem$Karyotype.band,
         min)

largest_order_kb_OS1 <-
  tapply(chr17_OS1_copy_number_outrem$order_17,
         chr17_OS1_copy_number_outrem$Karyotype.band,
         max)

positions_kb_OS1 <- data.frame(
  kb = names(smallest_order_kb_OS1),
  smallest_order = smallest_order_kb_OS1,
  largest_order = largest_order_kb_OS1
)

positions_kb_OS1$smallest_order <- positions_kb_OS1$smallest_order - 1

# centromere
centromere_row <- data.frame(
  kb = "",
  smallest_order = 323,
  largest_order = 376
)

positions_kb_OS1 <- rbind(positions_kb_OS1, centromere_row)

positions_kb_OS1$smallest_order[positions_kb_OS1$kb == "p13.3"] <- 1
positions_kb_OS1$smallest_order <- as.numeric(positions_kb_OS1$smallest_order)
positions_kb_OS1$largest_order <- as.numeric(positions_kb_OS1$largest_order)

positions_kb_OS1$mid <-
  (positions_kb_OS1$smallest_order + positions_kb_OS1$largest_order) / 2

smallest_order_kb_OS1 <- tapply(chr17_OS1_copy_number_outrem$order_17, chr17_OS1_copy_number_outrem$Karyotype.band, min)
largest_order_kb_OS1 <- tapply(chr17_OS1_copy_number_outrem$order_17, chr17_OS1_copy_number_outrem$Karyotype.band, max)

positions_kb_OS1 <- data.frame(
  kb = names(smallest_order_kb_OS1),
  smallest_order = smallest_order_kb_OS1,
  largest_order = largest_order_kb_OS1
)
positions_kb_OS1$smallest_order<- positions_kb_OS1$smallest_order -1 

#Manually insert a row for the centromere
centromere_row<-data.frame(kb = "", smallest_order = "323", largest_order = "376")
positions_kb_OS1<-rbind(positions_kb_OS1, centromere_row)

positions_kb_OS1$smallest_order[positions_kb_OS1$kb == "p13.3"] <- 1
positions_kb_OS1$smallest_order<-as.numeric(positions_kb_OS1$smallest_order)
positions_kb_OS1$largest_order<-as.numeric(positions_kb_OS1$largest_order)

#obtain label positions
positions_kb_OS1$mid<- (positions_kb_OS1$smallest_order + positions_kb_OS1$largest_order) / 2

light_kb<-c("p11.2", "p13.1", "p13.3", "q11.2", "q21.1", "q21.31", "q21.33", "q23.1", "q23.3", "q24.2" ,"q25.1" ,"q25.3" )
dark_kb<-c("p12" , "p13.2"  ,  "q12" , "q21.2" ,  "q21.32", "q22", "q23.2" ,  "q24.1" ,  "q24.3" , "q25.2") 

light_positions_kb_OS1<-positions_kb_OS1[positions_kb_OS1$kb %in% light_kb,]
dark_positions_kb_OS1<-positions_kb_OS1[positions_kb_OS1$kb %in% dark_kb,]
centromere_17OS1<-positions_kb_OS1[positions_kb_OS1$kb == "",]


light_kb_rectangles<-data.frame(
  xmin_chr17 = light_positions_kb_OS1$smallest_order,
  xmax_chr17 = light_positions_kb_OS1$largest_order,
  ymin_chr17 = -1.5,
  ymax_chr17 = -0.5
)


dark_kb_rectangles<-data.frame(
  xmin_chr17 = dark_positions_kb_OS1$smallest_order,
  xmax_chr17 = dark_positions_kb_OS1$largest_order,
  ymin_chr17 = -1.5,
  ymax_chr17 = -0.5
)


centromere_rectangles<-data.frame(
  xmin_chr17 = centromere_17OS1$smallest_order,
  xmax_chr17 = centromere_17OS1$largest_order,
  ymin_chr17 = -1.5,
  ymax_chr17 = -0.5
)



#Make plots
average_data_OS1_chr17 <- chr17_OS1_copy_number_outrem %>%
  group_by(order_17) %>%
  summarise(average_cn = mean(cn), .groups = 'drop')


chr17_OS1_copy_number_outrem <- chr17_OS1_copy_number_outrem %>%
  mutate(color = case_when(
    cn > 2 ~ "red",      # Red for cn > 2
    cn < 2 ~ "blue",     # Blue for cn < 2
    cn == 2 ~ "darkgrey"     # Grey for cn = 2
  )) 


darklables<-dark_positions_kb_OS1[dark_positions_kb_OS1$kb %in% c("p13.2", "p12", "q12", "q21.2", "q21.32", "q22"),]
lightlables<-light_positions_kb_OS1[light_positions_kb_OS1$kb %in% c("p11.2", "p13.1", "p13.3", "q11.2", "q21.31", "q21.33", "q23.3","q25.1" ,"q25.3" ) ,]



# Create the dot plot with jitter and average segments

OS1_cn_17_ploidy <- ggplot() + 
  # Add the jittered points
  geom_point(data = chr17_OS1_copy_number_outrem, 
             aes(x = order_17, y = cn, color = color), 
             size = 0.1, 
             position = position_jitter(width = 0.5, height = 0.5), 
             alpha = 1) +  
  # Add horizontal segments for average values
  geom_segment(data = average_data_OS1_chr17, aes(x = order_17, xend = order_17, y = average_cn, yend = average_cn), color = "black", size = 0.5) + 
  # Add points for average values for better visibility
  geom_point(data = average_data_OS1_chr17, aes(x = order_17, y = average_cn), color = "black", size = 0.001, shape = 0.01, fill = "black") + 
  # Add rectangles
  geom_rect(data = dark_kb_rectangles, 
            aes(xmin = xmin_chr17, xmax = xmax_chr17, ymin = ymin_chr17, ymax = ymax_chr17), 
            fill = "black", alpha = 1) +  
  geom_rect(data = light_kb_rectangles, 
            aes(xmin = xmin_chr17, xmax = xmax_chr17, ymin = ymin_chr17, ymax = ymax_chr17), 
            fill = "lightgrey", alpha = 1) + 
  geom_rect(data =  centromere_rectangles, 
            aes(xmin = xmin_chr17, xmax = xmax_chr17, ymin = ymin_chr17, ymax = ymax_chr17), 
            fill = "red", alpha = 1) + 
  geom_text(data = darklables, 
            aes(x = mid, y = -0.95, label = kb),  # Adjust y position as needed
            color = "white", size = 2.5, vjust = -0.5,fontface = "bold") +
  geom_text(data = lightlables, 
            aes(x = mid, y = -0.95, label = kb),  # Adjust y position as needed
            color = "white", size = 2.5, vjust = -0.5,fontface = "bold") +
  theme_classic() + 
  labs(x = "", y = "Copy number", title = "OS1 chromosome 17 copy number") + 
  theme(axis.text.x = element_blank(),           # Remove x-axis text
        axis.ticks.x = element_blank(),          # Remove x-axis ticks
        axis.line.x = element_blank(),            # Remove x-axis line
        axis.text.y = element_text(size = 12),   # Keep y-axis text
        plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks = 0:20) + 
  scale_x_continuous(expand = c(0, 0)) +  # Adjusted for numeric x-axis
  scale_color_identity()
print(OS1_cn_17_ploidy)
ggsave("/Users/catincaapostol/Desktop/Master's/wgs os1 analysis plots/OS1_cn_17_ploidy.jpg", plot = OS1_cn_17_ploidy, width = 12, height = 4, dpi = 300)


