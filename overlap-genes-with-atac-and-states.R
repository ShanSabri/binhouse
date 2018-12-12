### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/overlap-genes-with-atak-and-state/")



### LIBS
pacman::p_load(
  pheatmap, ggplot2,
  biomaRt, GenomicRanges,
  rtracklayer, tidyverse
)



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

feature_overlap <- function(peak, tss, extend_upstream = 20000, extend_downstream = 20000) {
  names(peak) <- seq(1:length(peak))
  window <- seq(1:length(tss))
  tss_ext <- extend(tss[window], extend_upstream, extend_downstream)
  names(tss_ext) <- window

  hits <- findOverlaps(tss_ext, peak, ignore.strand = FALSE)

  overlaps <- data.frame(
    gene_id = tss$gene_id[hits@from],
    gene_chr = chrom(tss)[hits@from],
    gene_tss = start(tss)[hits@from],
    gene_start = start(tss_ext)[hits@from],
    gene_end = end(tss_ext)[hits@from],
    peak_id = names(peak)[hits@to],
    peak_start = start(peak)[hits@to],
    peak_end = end(peak)[hits@to]
  )

  return(overlaps)
}

run_homer <- function(file, out, fa = " ~/Dropbox/Binhouse/mm9.fa", bg = NULL, dry = TRUE) {
  Sys.setenv(PATH = paste0(
    "/usr/local/opt/libxml2/bin:/Users/shansabri/miniconda3/bin:/Library/Frameworks/Python.framework/Versions/3.4/bin:",
    "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/ncbi/blast/bin:/Library/TeX/texbin:/usr/local/MacGPG2/bin:",
    "/opt/X11/bin:/Users/shansabri/bin:/Users/shansabri/edirect:/Users/shansabri/homer/bin/:/Users/shansabri/bin:/Users/shansabri/edirect"
  ))
  cmd <- sprintf("findMotifsGenome.pl %s %s %s -p 6 -size given", file, fa, out)
  if (!is.null(bg)) {
    cmd <- sprintf("findMotifsGenome.pl %s %s %s -bg %s -p 6 -size given", file, fa, out, bg)
  }

  print(cmd)

  if (dry == TRUE) {
    return(cmd)
  } else {
    try(system(cmd))
    return(TRUE)
  }
}




### ATAC DATA
dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/kostas_atac/"
beds <- list.files(dir, pattern = "*.bed", full.names = TRUE)
atac <- lapply(beds, bed_to_granges)
names(atac) <- basename(beds)
atac <- do.call(rbind.data.frame, lapply(atac, as.data.frame))
atac <- with(atac, GRanges(seqnames, IRanges(start, end))) # 402503 ATAC peaks




### CHROMATIN STATE DATA
state_annot <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/chromatin-states/STACKED-35-STATE-MODEL.rds")




### PROGRAM DATA
programs <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/Drop10.Aggregated.Programs.Gene.Master.List.rds")
programs_feat <- generate_gene_annotation_biomart(as.vector(unique(programs$Genes)))
programs_feat <- reduce_to_tss(programs_feat)
programs <- merge(programs, as.data.frame(programs_feat), by.x = "Genes", by.y = "symbol", all.y = TRUE)
programs <- split(programs, f = programs$Set)
programs <- lapply(programs, function(x) with(x, GRanges(seqnames, IRanges(start, end), gene_id = x$Genes)))




### COMPUTE OVERLAP OF PEAK GIVEN PROGRAM GENES
overlaps <- lapply(programs, function(p) feature_overlap(peak = atac, tss = p, extend_upstream = 20000, extend_downstream = 20000))
# lapply(overlaps, nrow)




### HOMER
homer <- lapply(seq_along(overlaps), function(o) {
  dir <- names(overlaps)[o]
  dir.create(dir)
  b <- with(overlaps[[o]], GRanges(gene_chr, IRanges(peak_start, peak_end)))
  f <- paste(dir, paste(dir, "bed", sep = "."), sep = "/")
  write.table(data.frame(seqnames = seqnames(b), starts = start(b) - 1, ends = end(b)),
    file = f, quote = FALSE, sep = "\t", row.names = F, col.names = F
  )
  run_homer(file = f, out = dir, bg = NULL, dry = FALSE)
})
