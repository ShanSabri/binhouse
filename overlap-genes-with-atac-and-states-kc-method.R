### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/overlap-genes-with-atac-and-state-kc-method/")



### LIBS
pacman::p_load(
  pheatmap, ggplot2, gplots,
  biomaRt, GenomicRanges,
  rtracklayer, tidyverse, bedr, marge
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

state_overlap <- function(peak, states, extend_upstream = 0, extend_downstream = 0) {
  names(peak) <- seq(1:length(peak))
  window <- seq(1:length(states))
  states_ext <- extend(states[window], extend_upstream, extend_downstream)
  names(states_ext) <- window
  
  hits <- findOverlaps(states_ext, peak, ignore.strand = FALSE)
  
  overlaps <- data.frame(
    state = states$state[hits@from],
    state_chr = chrom(states)[hits@from],
    state_start = start(states_ext)[hits@from],
    state_end = end(states_ext)[hits@from],
    peak_id = names(peak)[hits@to],
    peak_start = start(peak)[hits@to],
    peak_end = end(peak)[hits@to]
  )
  
  overlaps <- with(overlaps, GRanges(state_chr, IRanges(peak_start, peak_end)))
  
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




# ### ATAC DATA
# dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/kostas_atac/"
# beds <- list.files(dir, pattern = "*.bed", full.names = TRUE)
# atac <- lapply(beds, bed_to_granges)
# names(atac) <- basename(beds)
# 
# atac$MEF_ATACseq.peaks.bed <- NULL
# atac$`48h_ATACseq.peaks.bed` <- NULL
# atac$preiPSC1_ATACseq.peaks.bed <- NULL
# atac$ESC_ATACseq.peaks.bed <- NULL
# 
# 
# atac <- do.call(rbind.data.frame, lapply(atac, as.data.frame))
# atac <- with(atac, GRanges(seqnames, IRanges(start, end))) # 402503 ATAC peaks
# atac_bed <- data.frame(seqnames=seqnames(atac), starts=start(atac)-1, ends=end(atac))
# # write.table(atac_bed, file=paste(getwd(), "ALL-MEF-ATAC-PEAKS.bed", sep = "/"), quote=F, sep="\t", row.names=F, col.names=F)
# # write.table(atac_bed, file=paste(getwd(), "ALL-48H-PIPS-ATAC-PEAKS.bed", sep = "/"), quote=F, sep="\t", row.names=F, col.names=F)
# write.table(atac_bed, file=paste(getwd(), "ALL-ESC-ATAC-PEAKS.bed", sep = "/"), quote=F, sep="\t", row.names=F, col.names=F)
# 
# 
# 
# 
# ### CHROMATIN STATE DATA
# state_annot <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/chromatin-states/STACKED-35-STATE-MODEL.rds")
# 
# 
# 
# 
# ### LOOK ONLY AT PEAKS THAT OVERLAP CHROMATIN STATES
# # int <- state_annot[ state_annot$state %in% paste0("E", seq(5, 10, 1)) ]
# # int <- state_annot[ state_annot$state %in% paste0("E", seq(11, 12, 1)) ]
# int <- state_annot[ state_annot$state %in% paste0("E", seq(13, 18, 1)) ]
# atac <- state_overlap(peak = atac, states = int, extend_upstream = 0, extend_downstream = 0)
# 
# 
# ### EXTEND PEAKS OUT AND MERGE
# summary(width(atac))
# ggplot(data.frame(w = width(atac)), aes(w)) + geom_density()
# x = as.data.frame(atac)
# x$start[x$width == 2] <- x$start[x$width == 2] - 200
# x$end[x$width == 2] <- x$end[x$width == 2] + 200
# x$width <- x$end - x$start
# 
# options(scipen=999)
# merged_peaks <- bedr.merge.region(
#   x = paste0(x$seqnames, ":", x$start, "-", x$end),
#   distance = 2,
#   list.names = TRUE,
#   number = FALSE,
#   stratify.by = NULL,
#   check.zero.based = TRUE,
#   check.chr = TRUE,
#   check.valid = TRUE,
#   check.sort = TRUE,
#   verbose = TRUE
# ) #  Collapsing 237127 --> 168541 regions
# bed <- convert2bed(merged_peaks)
# # write.table(bed, file=paste(getwd(), "ALL-MEF-ATAC-PEAKS-OVERLAP-ENH-5-10.bed", sep = "/"), quote=F, sep="\t", row.names=F, col.names=F)
# # write.table(bed, file=paste(getwd(), "ALL-48H-PIPS-ATAC-PEAKS-OVERLAP-ENH-11-12.bed", sep = "/"), quote=F, sep="\t", row.names=F, col.names=F)
# write.table(bed, file=paste(getwd(), "ALL-ESC-ATAC-PEAKS-OVERLAP-ENH-13-18.bed", sep = "/"), quote=F, sep="\t", row.names=F, col.names=F)
# 

####################################################################################################
## STG 2 -- RUN PROGRAM OVERLAP ONCE WE HAVE THE BED FILES




### PROGRAM DATA
programs <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/Drop10.Aggregated.Programs.Gene.Master.List.rds")
programs <- droplevels.data.frame(subset(programs, Set %in% c("CDH1_1"))) # subset) , "CDH1_8", "TC1_4.02"
programs_feat <- generate_gene_annotation_biomart(as.vector(unique(programs$Genes)))
programs_feat <- reduce_to_tss(programs_feat)
programs <- merge(programs, as.data.frame(programs_feat), by.x = "Genes", by.y = "symbol", all.y = TRUE)
programs <- split(programs, f = programs$Set)
programs <- lapply(programs, function(x) with(x, GRanges(seqnames, IRanges(start, end), gene_id = x$Genes)))




### READ IN BED AND CONVERT TO GRANGES 
bed <- read.table(file=paste(getwd(), "ALL-ESC-ATAC-PEAKS-OVERLAP-ENH-13-18.bed", sep = "/"), sep="\t")
atac <- GRanges(seqnames = bed$V1, ranges = IRanges(start = bed$V2, end = bed$V3))




### COMPUTE OVERLAP OF PEAK GIVEN PROGRAM GENES
overlaps <- lapply(programs, function(p) feature_overlap(peak = atac, tss = p, extend_upstream = 20000, extend_downstream = 20000))




### HOMER
ID <- "ESC"
homer <- lapply(seq_along(overlaps), function(o) {
  
  o = 1
  
  dir <- paste(names(overlaps)[o], ID, sep = "_")
  dir.create(dir, showWarnings = FALSE)
  b <- with(overlaps[[o]], GRanges(gene_chr, IRanges(peak_start, peak_end)))
  df <- data.frame(seqnames = seqnames(b), starts = start(b) - 1, ends = end(b))
  f <- paste(dir, paste(dir, "bed", sep = "."), sep = "/")
  write.table(df, file = f, quote = FALSE, sep = "\t", row.names = F, col.names = F)
  run_homer(file = f, out = dir, dry = TRUE)
})
getwd()



# 
# 
# ### VISUALIZE AGGREGRATED MOTIF ENRICHMENT
# dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/overlap-genes-with-atac-and-state-w-bg/"
# res <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
# all_motif_results <- do.call(rbind.data.frame, lapply(seq_along(res), function(x) {
#   # denovo <- read_denovo_results(path = x, homer_dir = TRUE)
#   known <- read_known_results(path = res[[x]], homer_dir = TRUE)
#   known$ID <- basename(res)[x]
#   known <- subset(known, database == "Homer")
#   known
# }))
# saveRDS(all_motif_results, compress = TRUE, paste(dir, "RESULTS.rds", sep = "/"))
# 
# # WHAT ARE THE TOP 10 MOTIFS FOR EACH GENE SET?
# all_motif_results %>%
#   group_by(ID) %>%
#   top_n(log_p_value, n = 25) %>%
#   dplyr::select(., motif_name, motif_family, experiment, accession, log_p_value) %>%
#   ungroup() %>%
#   write_tsv("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/overlap-genes-with-atac-and-state-w-bg/TOP25-MOTIFS-BY-PROGRAM.txt",
#             na = "NA",  quote_escape = "double")
# 
# 
# motif_df <- reshape2::dcast(all_motif_results, motif_name + motif_family + accession + consensus ~ ID,
#                             value.var = "log_p_value", fun.aggregate = max
# )
# row.names(motif_df) <- paste(motif_df$motif_name, motif_df$motif_family, motif_df$accession, motif_df$consensus, sep = "/")
# motif_df$motif_name <- NULL
# motif_df$motif_family <- NULL
# motif_df$accession <- NULL
# motif_df$consensus <- NULL
# 
# programs <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/Drop10.Aggregated.Programs.Gene.Master.List.rds")
# order_of_progs <- unique(programs[, c("Order", "Set")])
# motif_df <- motif_df[, as.vector(order_of_progs$Set)]
# motif_df$max_p_val <- apply(motif_df, 1, max)
# 
# maxs <- seq(1, 20, by = 2)
# pdf(paste(getwd(), paste0("ENRICHED_MOTIFS_MAX_THRESHOLD.pdf"), sep = "/"), height = 15, width = 10)
# lapply(maxs, function(x) {
#   tmp <- subset(motif_df, max_p_val >= x)
#   tmp$max_p_val <- NULL
#   # heatmap.2(as.matrix(log10(tmp + 1)),
#   heatmap.2(as.matrix(tmp),
#             Rowv = TRUE,
#             Colv = FALSE,
#             trace = "none",
#             na.color = "grey",
#             margins = c(8, 20),
#             col = cm.colors(255),
#             main = paste("Max log_p_value >=", x)
#   )
# })
# dev.off()
