
maf_file <- "/Users/catincaapostol/Desktop/Master's/project/wgs os1/NRHOS1_SNVcalling_filtered.maf.gz"

lines <- readLines(gzfile(maf_file))
header_idx <- which(grepl("Hugo_Symbol", lines))[1]
if (is.na(header_idx)) stop("Header line with 'Hugo_Symbol' not found in the file.")

header_line <- lines[header_idx]

col_names <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
library(data.table)
install.packages('R.utils')
library(R.utils)
snv_calling_df <- fread(
  maf_file,
  sep = "\t",
  skip = "Hugo_Symbol",
  header = TRUE,
  colClasses = "character",
  quote = "",
  fill = TRUE,
  showProgress = FALSE
)

cat("Rows (variants):", nrow(snv_calling_df), "\n")
cat("Columns:", ncol(snv_calling_df), "\n")
stopifnot(ncol(snv_calling_df) == length(col_names))
sort(unique(snv_calling_df$Chromosome))

#😘😀😆😄👏🏻🤝👎🏻❤️🧡💛💚🩵💙💜🩷🖤🩶🤍🤎
snv_prep_df <- snv_calling_df[!grepl(",", snv_calling_df$TLOD, fixed = TRUE), ]
snv_filter_ready_df <- snv_prep_df

snv_filter_ready_df$MMQ <- as.numeric(sub(".*,(.*)\\]", "\\1", snv_filter_ready_df$MMQ))
snv_filter_ready_df$MBQ <- as.numeric(sub(".*,(.*)\\]", "\\1", snv_filter_ready_df$MBQ))
#convert numeric columns 
snv_filter_ready_df$DP <- as.numeric(snv_filter_ready_df$DP)
snv_filter_ready_df$TLOD <- as.numeric(snv_filter_ready_df$TLOD)
snv_filter_ready_df$GERMQ <- as.numeric(snv_filter_ready_df$GERMQ)
snv_filter_ready_df$t_alt_count <- as.numeric(snv_filter_ready_df$t_alt_count)

  # if DP>=20 and t_alt_count>=5 
    # if TLOD >=8
      #if MMQ>=30 
        #if MBQ>=25
          #if GERMQ >=30 
            #keep row with all columns to a new dataframe called snv_priority_calls 

snv_priority_calls1 <- snv_filter_ready_df[
  snv_filter_ready_df$DP >= 17, ]


snv2_priority_calls2 <- snv_priority_calls1[
  snv_priority_calls1$t_alt_count >= 5, 
]
rm(snv_priority_calls1)
snv3_priority_calls3 <- snv2_priority_calls2[
  snv2_priority_calls2$TLOD >= 6, 
]
rm(snv2_priority_calls2)
snv4_priority_calls4 <- snv3_priority_calls3[
  snv3_priority_calls3$MMQ >= 20, 
]
rm(snv3_priority_calls3)
snv5_priority_calls5 <- snv4_priority_calls4[
  snv4_priority_calls4$MBQ >= 20, 
]
rm(snv4_priority_calls4)
snv6_priority_calls6 <- snv5_priority_calls5[
  snv5_priority_calls5$GERMQ >= 20, 
]
rm(snv5_priority_calls5)
snv_final_df <- snv6_priority_calls6[
  snv6_priority_calls6$Variant_Classification %in% c(
    "Missense_Mutation",
    "In_Frame_Ins",
    "Splice_Site",
    "Nonsense_Mutation",
    "Frame_Shift_Ins",
    "Frame_Shift_Del",
    "In_Frame_Del",
    "DE_NOVO_START_IN_FRAME",
    "DE_NOVO_START_OUT_FRAME",
    "Nonstop_Mutation",
    "START_CODON_SNP",
    "START_CODON_INS"
  ),
]
rownames(snv_final_df) <- NULL
rm(snv6_priority_calls6)
snv_final_df <- snv_final_df[
  snv_final_df$tumor_f >= 0., 
]
#snv1_mid_confidence <- snv_filter_ready_df[
  #snv_filter_ready_df$DP >= 10 & snv_filter_ready_df$DP <= 20, ]

#snv2_mid_confidence <- snv1_mid_confidence[
#  snv1_mid_confidencef$t_alt_count >= 5, ]
#rm(snv1_mid_confidence)
#snv3_mid_confidence <- snv2_mid_confidence[
#  snv2_mid_confidencef$TLOD >= 5 & snv2_mid_confidencef$TLOD <= 8]
#rm(snv2_mid_confidence)
#snv4_mid_confidence <- snv3_mid_confidence[
#  snv3_mid_confidencef$TLOD >= 5 & snv3_mid_confidencef$TLOD <= 8]
#rm(snv3_mid_confidence)
