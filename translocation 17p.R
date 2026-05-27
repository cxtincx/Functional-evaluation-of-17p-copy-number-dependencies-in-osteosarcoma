library (ggplot2)
translocation_vcf <- "/Users/catincaapostol/Desktop/Master's/project/wgs os1/candidateSV.vcf.gz"
cran_pkgs <- c("data.table", "circlize")
bioc_pkgs  <- c("VariantAnnotation", "GenomicRanges", "GenomicFeatures",
                "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db", "AnnotationDbi")
install_if_missing <- function(pkgs, bioconductor=FALSE) {
  if(length(pkgs)==0) return()
  for(pkg in pkgs) {
    if(!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) {
      if(bioconductor) {
        if(!suppressWarnings(requireNamespace("BiocManager", quietly = TRUE))) {
          install.packages("BiocManager", repos="https://cloud.r-project.org")
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, repos="https://cloud.r-project.org")
      }
    }
  }
}

install_if_missing(cran_pkgs, bioconductor = FALSE)
install_if_missing(bioc_pkgs, bioconductor = TRUE)

# load libraries
library(data.table)
library(circlize)
library(VariantAnnotation)
library(GenomicRanges)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(AnnotationDbi)

norm_chr <- function(x) {
  x <- as.character(x)
  x <- sub("^chr", "", x, ignore.case = TRUE)
  paste0("chr", x)
}

parse_alt_chrpos <- function(altstr) {
  if(is.na(altstr) || altstr == "") return(NULL)
  # find first occurrence of pattern like chrX:12345 (allow underscores/letters)
  m <- regmatches(altstr, regexpr("([A-Za-z0-9_\\.\\-]+):([0-9]+)", altstr))
  if(length(m)==0 || m=="") return(NULL)
  parts <- strsplit(m, ":")[[1]]
  list(chr = parts[1], pos = as.integer(parts[2]))
}


get_info_field <- function(info_df, field_name, n) {
  if(!(field_name %in% colnames(info_df))) return(rep(NA_character_, n))
  col <- info_df[[field_name]]
  if(is(col, "List") || is.list(col) || inherits(col, "CompressedIntegerList") || inherits(col, "IntegerList") || inherits(col, "CharacterList")) {
    res <- sapply(col, function(x) { if(length(x)==0) NA_character_ else as.character(x[1]) }, USE.NAMES = FALSE)
    if(length(res) != n) res <- rep(NA_character_, n)
    return(res)
  } else {
    res <- as.character(col)
    if(length(res) != n) res <- rep(NA_character_, n)
    return(res)
  }
}


map_entrez_to_symbol <- function(entrez_ids) {
  if(length(entrez_ids)==0) return(character(0))
  syms <- mapIds(org.Hs.eg.db, keys = as.character(entrez_ids),
                 column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  return(syms)
}


canonical_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")


message("Reading VCF (may take a moment)...")
vcf <- readVcf(translocation_vcf, "hg38")
rr <- rowRanges(vcf)
nrec <- length(rr)
message("Total VCF records: ", nrec)

info_df <- info(vcf)
alts_list <- alt(vcf)



alts <- sapply(alts_list, function(x) {
  if(length(x)==0) return(NA_character_)
  paste0(as.character(x), collapse=",")
}, USE.NAMES = FALSE)



END_field    <- get_info_field(info_df, "END", nrec)
CHR2_field   <- get_info_field(info_df, "CHR2", nrec)
MATEID_field <- get_info_field(info_df, "MATEID", nrec)
SVTYPE_field <- get_info_field(info_df, "SVTYPE", nrec)
TYPE_field   <- get_info_field(info_df, "TYPE", nrec)

svtype_vec <- SVTYPE_field
svtype_vec[is.na(svtype_vec) | svtype_vec == "NA"] <- TYPE_field[is.na(svtype_vec) | svtype_vec == "NA"]
svtype_vec[is.na(svtype_vec)] <- NA_character_



chr1 <- as.character(seqnames(rr))
pos1 <- start(rr)
idv  <- names(rr)
filter_vec <- as.character(fixed(vcf)$FILTER)
if(length(filter_vec) != nrec) filter_vec <- rep(NA_character_, nrec)



chr2_parsed <- rep(NA_character_, nrec)
pos2_parsed <- rep(NA_integer_, nrec)



bnd_idx <- which(grepl("\\[|\\]", alts))
if(length(bnd_idx) > 0) {
  for(i in bnd_idx) {
    parsed <- parse_alt_chrpos(alts[i])
    if(!is.null(parsed)) {
      chr2_parsed[i] <- parsed$chr
      pos2_parsed[i] <- parsed$pos
    }
  }
}



has_chr2 <- which(!is.na(CHR2_field) & CHR2_field != "")
if(length(has_chr2) > 0) {
  for(i in has_chr2) {
    if(is.na(chr2_parsed[i])) {
      chr2_parsed[i] <- CHR2_field[i]
      # if END exists use it for pos
      if(!is.na(END_field[i]) && END_field[i] != "") pos2_parsed[i] <- as.integer(END_field[i])
    }
  }
}



if(any(!is.na(MATEID_field))) {
  id_to_idx <- setNames(seq_len(nrec), idv)
  for(i in seq_len(nrec)) {
    if(is.na(chr2_parsed[i]) && !is.na(MATEID_field[i]) && MATEID_field[i] != "") {
      mids <- strsplit(MATEID_field[i], ",")[[1]]
      mate_first <- mids[1]
      if(mate_first %in% names(id_to_idx)) {
        j <- id_to_idx[[mate_first]]
        chr2_parsed[i] <- as.character(seqnames(rr)[j])
        pos2_parsed[i] <- start(rr)[j]
      }
    }
  }
}



for(i in seq_len(nrec)) {
  if(is.na(chr2_parsed[i]) && !is.na(END_field[i]) && END_field[i] != "") {
    chr2_parsed[i] <- chr1[i]
    pos2_parsed[i] <- as.integer(END_field[i])
  }
}

qual_vec <- fixed(vcf)$QUAL

sv_dt <- data.table(
  id = ifelse(is.null(idv), paste0("rec", seq_len(nrec)), idv),
  chr1 = norm_chr(chr1),
  pos1 = as.integer(pos1),
  chr2 = ifelse(is.na(chr2_parsed), NA_character_, norm_chr(chr2_parsed)),
  pos2 = as.integer(pos2_parsed),
  svtype = svtype_vec,
  ALT = alts,
  stringsAsFactors = FALSE
)



message("Parsed partner coords for ", sum(!is.na(sv_dt$chr2) & !is.na(sv_dt$pos2)), " / ", nrec, " records.")


sv_dt[, chr1_short := sub("^chr", "", chr1, ignore.case = TRUE)]
sv_dt[, chr2_short := ifelse(is.na(chr2), NA_character_, sub("^chr", "", chr2, ignore.case = TRUE))]

trans_candidates <- sv_dt[!is.na(chr2) & !is.na(pos2) & chr1_short != chr2_short]

message("Inter-chromosomal candidates (all contigs): ", nrow(trans_candidates))


trans_canonical <- trans_candidates[chr1 %in% canonical_chrs & chr2 %in% canonical_chrs]

message("Inter-chromosomal candidates on canonical chromosomes (chr1-22,X,Y,M): ", nrow(trans_canonical))
if(nrow(trans_canonical) == 0) {
  message("No canonical-chromosome inter-chromosomal events found. There are candidates on non-canonical contigs (decoy/random).")
  fwrite(trans_candidates, file = "manta_translocation_candidates_allcontigs.tsv", sep = "\t")
  stop("No canonical translocations to plot. See manta_translocation_candidates_allcontigs.tsv for details.")
}


txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)  


bp1_gr <- GRanges(seqnames = trans_canonical$chr1, ranges = IRanges(start = trans_canonical$pos1, width = 1))
bp2_gr <- GRanges(seqnames = trans_canonical$chr2, ranges = IRanges(start = trans_canonical$pos2, width = 1))

ov1 <- findOverlaps(bp1_gr, genes_gr)
ov2 <- findOverlaps(bp2_gr, genes_gr)


genes1 <- rep(NA_character_, length(bp1_gr))
genes2 <- rep(NA_character_, length(bp2_gr))

if(length(ov1) > 0) {
  hits1 <- split(subjectHits(ov1), queryHits(ov1))
  for(q in names(hits1)) {
    entrez_ids <- names(genes_gr)[hits1[[q]]]  
    syms <- map_entrez_to_symbol(entrez_ids)
    genes1[as.integer(q)] <- paste(unique(na.omit(syms)), collapse = ";")
  }
}
if(length(ov2) > 0) {
  hits2 <- split(subjectHits(ov2), queryHits(ov2))
  for(q in names(hits2)) {
    entrez_ids <- names(genes_gr)[hits2[[q]]]
    syms <- map_entrez_to_symbol(entrez_ids)
    genes2[as.integer(q)] <- paste(unique(na.omit(syms)), collapse = ";")
  }
}


nearest_window <- 10000L
no1 <- which(is.na(genes1) | genes1=="")
no2 <- which(is.na(genes2) | genes2=="")

if(length(no1) > 0) {
  near1 <- nearest(bp1_gr[no1], genes_gr)
  for(i in seq_along(no1)) {
    nid <- near1[i]
    if(!is.na(nid)) {
      distv <- distance(bp1_gr[no1[i]], genes_gr[nid])
      if(is.na(distv) || distv <= nearest_window) {
        eid <- names(genes_gr)[nid]
        genes1[no1[i]] <- map_entrez_to_symbol(eid)
      }
    }
  }
}
if(length(no2) > 0) {
  near2 <- nearest(bp2_gr[no2], genes_gr)
  for(i in seq_along(no2)) {
    nid <- near2[i]
    if(!is.na(nid)) {
      distv <- distance(bp2_gr[no2[i]], genes_gr[nid])
      if(is.na(distv) || distv <= nearest_window) {
        eid <- names(genes_gr)[nid]
        genes2[no2[i]] <- map_entrez_to_symbol(eid)
      }
    }
  }
}


out_dt <- data.table(
  id = trans_canonical$id,
  chr1 = trans_canonical$chr1, pos1 = trans_canonical$pos1, genes1 = ifelse(is.na(genes1),"",genes1),
  chr2 = trans_canonical$chr2, pos2 = trans_canonical$pos2, genes2 = ifelse(is.na(genes2),"",genes2),
  svtype = trans_canonical$svtype,
  FILTER = trans_canonical$FILTER,
  ALT = trans_canonical$ALT
)
fwrite(out_dt, file = "manta_translocations_annotated.tsv", sep = "\t")
message("Wrote annotated translocations to manta_translocations_annotated.tsv (rows = ", nrow(out_dt), ").")


# ===== Full replacement block - zoom in to 17p 
if(!requireNamespace("circlize", quietly = TRUE)) install.packages("circlize", repos = "https://cloud.r-project.org")
if(!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
if(!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer", repos = "https://cloud.r-project.org")

library(circlize)
library(data.table)
library(grid)
library(RColorBrewer)


if(!exists("out_dt")) {
  if(file.exists("manta_translocations_annotated.tsv")) {
    out_dt <- fread("manta_translocations_annotated.tsv", data.table = FALSE)
  } else stop("out_dt not found and manta_translocations_annotated.tsv not present. Run the parser first.")
}


out_dt$chr1 <- as.character(out_dt$chr1)
out_dt$chr2 <- as.character(out_dt$chr2)
out_dt$pos1 <- as.integer(as.numeric(out_dt$pos1))
out_dt$pos2 <- as.integer(as.numeric(out_dt$pos2))


canonical_chrs <- c(paste0("chr", 1:22), "chrX", "chrY", "chrM")


valid_idx <- which(!is.na(out_dt$chr1) & !is.na(out_dt$chr2) &
                     out_dt$chr1 %in% canonical_chrs & out_dt$chr2 %in% canonical_chrs &
                     !is.na(out_dt$pos1) & !is.na(out_dt$pos2))
plot_all <- out_dt[valid_idx, , drop = FALSE]
if(nrow(plot_all) == 0) stop("No canonical translocation events found to plot.")

focal <- "chr17"
focal_start <- 14000000
focal_end   <- 22000000

sel_idx <- which(
  (plot_all$chr1 == focal & plot_all$pos1 >= focal_start & plot_all$pos1 <= focal_end) |
    (plot_all$chr2 == focal & plot_all$pos2 >= focal_start & plot_all$pos2 <= focal_end)
)
sel_dt  <- plot_all[sel_idx, , drop = FALSE]
if(nrow(sel_dt) == 0) stop("No translocations involving chr17 found in the parsed data.")


region1_all <- data.frame(chr = sel_dt$chr1, start = sel_dt$pos1, end = sel_dt$pos1, stringsAsFactors = FALSE)
region2_all <- data.frame(chr = sel_dt$chr2, start = sel_dt$pos2, end = sel_dt$pos2, stringsAsFactors = FALSE)


partner_chr <- ifelse(
  sel_dt$chr1 == focal & sel_dt$pos1 >= focal_start & sel_dt$pos1 <= focal_end,
  sel_dt$chr2,
  sel_dt$chr1
)
partners_unique <- sort(unique(partner_chr))
partners_unique <- partners_unique[partners_unique != focal]


nparts <- max(1, length(partners_unique))
if(nparts <= 8) {
  pal_base <- if("Dark2" %in% rownames(brewer.pal.info)) "Dark2" else "Set2"
  cols_palette <- brewer.pal(min(nparts, brewer.pal.info[pal_base, "maxcolors"]), pal_base)[1:nparts]
} else {
  cols_palette <- colorRampPalette(c("#006400", "#4B0082", "#2F4F4F", "#3B4994", "#2A6F97"))(nparts)
}
partner_colors <- setNames(cols_palette, partners_unique)


partner_hex <- sapply(partner_colors, function(x) {
  grDevices::rgb(t(grDevices::col2rgb(x))/255)
}, USE.NAMES = TRUE)


row_partner_colors_hex <- sapply(partner_chr, function(ch) {
  if(!is.na(ch) && ch %in% names(partner_hex)) partner_hex[ch] else grDevices::rgb(0.4,0.4,0.4)
}, USE.NAMES = FALSE)


draw_circos_page <- function(page_title, region1_local, region2_local, row_colors_hex, focal_chr) {
  used_chrs <- unique(c(focal_chr, region1_local$chr, region2_local$chr))
  used_chrs <- intersect(canonical_chrs, used_chrs)
  used_chrs <- c(focal_chr, setdiff(used_chrs, focal_chr))
  circos.clear()
  circos.par(start.degree = 90, gap.degree = 2, cell.padding = c(0,0,0,0), track.margin = c(0,0))
  circos.initializeWithIdeogram(species = "hg38", chromosome.index = used_chrs)
  circos.trackPlotRegion(track.height = 0.05, ylim = c(0,1), bg.border = NA,
                         panel.fun = function(region, value, ...) {
                           circos.text(mean(CELL_META$xlim), 0.5, CELL_META$sector.index, cex = 2, facing = "inside")
                         })
  
  
  nrows_local <- nrow(region1_local)
  for(i in seq_len(nrows_local)) {
    if(region1_local$chr[i] %in% used_chrs && region2_local$chr[i] %in% used_chrs) {
      circos.genomicLink(region1_local[i, , drop = FALSE], region2_local[i, , drop = FALSE],
                         col = row_colors_hex[i], border = NA)
    }
  }
  
  
  legend_labels <- names(partner_colors)
  legend_cols <- partner_colors
  if(length(legend_labels) > 0) {
    par(xpd = TRUE)
    legend("topright", legend = legend_labels, fill = legend_cols, border = NA, bty = "n", cex = 2, inset = c(0.01, 0.02))
  }
  
 
  grid::grid.text(label = page_title, x = 0.5, y = grid::unit(0.975, "npc"),
                  gp = grid::gpar(fontsize = 14, fontface = "bold"))
  
  circos.clear()
}


pdf("/Users/catincaapostol/Desktop/Master's/project/wgs os1/chr17p12-p11.2_translocations.pdf", width = 11, height = 10.5)
message("Writing chr17p12-p11.2_translocations.pdf")


d <- draw_circos_page(
  page_title = paste0("chr17p12-p11.2 translocations"),
  region1_local = region1_all,
  region2_local = region2_all,
  row_colors_hex = row_partner_colors_hex,
  focal_chr = focal
)
print(d)
dev.off()
message("Done — output file: chr17_translocations.pdf")

