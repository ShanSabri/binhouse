### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/") # set working directory



### LIBS
pacman::p_load(pheatmap, biomaRt, GenomicRanges, rtracklayer, tidyverse)



### FXNS
generate_gene_annotation_biomart <- function(g) {
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl", host = "may2012.archive.ensembl.org"))
  mapping <- getBM(
    filters = "mgi_symbol",
    attributes = c("chromosome_name", "start_position", "end_position", "strand", "mgi_symbol"),
    values = g, mart = mart
  )
  mapping <- subset(mapping, chromosome_name %in% c(seq(1, 19), "X", "Y"))
  gene <- GRanges(
    seqnames = Rle(paste0("chr", mapping$chromosome_name)),
    ranges = IRanges(start = mapping$start_position, end = mapping$end_position, names = seq(1, nrow(mapping))),
    strand = Rle(mapping$strand),
    symbol = toupper(mapping$mgi_symbol)
  )
  return(gene)
}

reduce_to_tss <- function(x) {
  return(resize(x, 1, fix = "start", use.names = TRUE, ignore.strand = FALSE))
}

bed_to_granges <- function(file) {
  df <- read.table(file, header = F, stringsAsFactors = F)
  if (length(df) > 6) df <- df[, -c(7:length(df))]
  if (length(df) < 3) stop("File has less than 3 columns")
  header <- c("chr", "start", "end", "id", "score", "strand")
  names(df) <- header[1:length(names(df))]

  if ("strand" %in% colnames(df)) {
    df$strand <- gsub(pattern = "[^+-]+", replacement = "*", x = df$strand)
  }

  if (length(df) == 3) {
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df) == 4) {
    gr <- with(df, GRanges(chr, IRanges(start, end), id = id))
  } else if (length(df) == 5) {
    gr <- with(df, GRanges(chr, IRanges(start, end), id = id, score = score))
  } else if (length(df) == 6) {
    gr <- with(df, GRanges(chr, IRanges(start, end), id = id, score = score, strand = strand))
  }
  return(gr)
}

extend <- function(x, upstream = 0, downstream = 0) {
  if (any(strand(x) == "*")) warning("'*' ranges were treated as '+'")
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  x
}

feature_overlap <- function(blocks, tss, extend_upstream = 20000, extend_downstream = 20000) {
  names(blocks) <- seq(1:length(blocks))
  window <- seq(1:length(tss))
  tss_ext <- extend(tss[window], extend_upstream, extend_downstream)
  names(tss_ext) <- window

  hits <- findOverlaps(tss_ext, blocks, ignore.strand = FALSE)

  overlaps <- data.frame(
    gene_id = tss$symbol[hits@from],
    gene_chr = chrom(tss)[hits@from],
    gene_tss = start(tss)[hits@from],
    gene_start = start(tss_ext)[hits@from],
    gene_end = end(tss_ext)[hits@from],
    peak_id = names(blocks)[hits@to],
    peak_start = start(blocks)[hits@to],
    peak_end = end(blocks)[hits@to]
  )

  return(overlaps)
}



### DATA
# - peak data
dir <- "~/dir/to/peaks"
beds <- list.files(dir, pattern = "*.bed", full.names = TRUE)
peaks <- lapply(beds, bed_to_granges)
names(peaks) <- basename(beds)

# - gene data
genes <- row.names(readRDS("~/dir/to/count_matrix"))
gene_feat <- generate_gene_annotation_biomart(genes)
gene_tss <- reduce_to_tss(gene_feat)



### COMPUTE OVERLAP GIVEN WINDOW
overlaps <- lapply(peaks, function(p) feature_overlap(blocks = p, tss = gene_tss, extend_upstream = 20000, extend_downstream = 20000))



### OUTPUT
out_list <- lapply(overlaps, function(x) {
  x %>%
    group_by(gene_id) %>%
    summarise(n_peaks = n()) %>%
    arrange(desc(n_peaks))
})
saveRDS(out_list, compress = TRUE, paste(getwd(), "GENES-BY-PEAKS-LIST.rds", sep = "/"))

out <- Reduce(function(x1, x2) merge(x1, x2, by = "gene_id", all = TRUE), out_list)
names(out) <- c("gene_id", names(out_list))
row.names(out) <- out$gene_id
out$gene_id <- NULL
out[is.na(out)] <- 0
saveRDS(out, compress = TRUE, paste(getwd(), "GENES-BY-PEAKS.rds", sep = "/"))


### VISUALIZE
library(pheatmap)
pheatmap(out,
  fontsize_col = 16,
  show_rownames = FALSE,
  filename = paste(getwd(), "HEATMAP-GENES-BY-PEAKS-RAW.png", sep = "/"),
  height = 20,
  width = 7
)

pheatmap(out[apply(out, 1, function(x) sum(x > 0)) < 5, ],
  fontsize_col = 16,
  fontsize_row = 3,
  filename = paste(getwd(), "HEATMAP-GENES-BY-PEAKS-FILT.png", sep = "/"),
  height = 25,
  width = 7
)
