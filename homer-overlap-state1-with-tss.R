################################################
## OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/homer/1111_3333/State1-overlap/"); getwd()
pacman::p_load(biomaRt, tidyverse, GenomicRanges)



################################################
## DEFINE PROGRAMS
gene_classes <- read.table("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/Clustering/gene-lists/DropAll.SharedGeneList.BinClass.and.PWM.txt", header = TRUE, sep = "\t")
gene_classes$PWMBin.Avg.Cut <- cut(gene_classes$PWMBin.Avg, breaks = 10, labels = FALSE)
gene_classes <- subset(gene_classes, Index %in% c(1111, 3333))
gene_classes$to_homer_group <- paste(gene_classes$Index, gene_classes$PWMBin.Avg.Cut, sep = "_")
tmp1 <- split(gene_classes, f = gene_classes$Index)
tmp2 <- split(gene_classes, f = gene_classes$to_homer_group); tmp2$`1111_3` <- NULL
programs <- c(tmp1, tmp2)




################################################
## GET PROGRAM INFO (TSS, CHR, ETC)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

get_gene_info <- function(x) {
  mapping <- getBM(filters = "mgi_symbol", 
                   attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "end_position", "start_position", "strand"),  
                   values = as.vector(x$Genes), 
                   mart = mart)
  message(paste("lookup loss percent", ((length(as.vector(x$Genes)) - nrow(mapping)) / length(as.vector(x$Genes))) * 100))

  to_return <- GenomicRanges::GRanges(
    seqnames = Rle(paste0("chr", mapping$chromosome_name)),
    ranges = IRanges(start = mapping$start_position, end = mapping$end_position, names = seq(1, nrow(mapping))),
    strand = Rle(mapping$strand),
    symbol = toupper(mapping$mgi_symbol)
  )
  
  return(to_return)
}

reduce_to_tss <- function(x) {
  return(resize(x, 1, fix = "start", use.names = TRUE, ignore.strand = FALSE))
}

extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*")){ warning("'*' ranges were treated as '+'")}
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  # trim(x)
  return(x)
}

run_homer <- function(file, fa, bg, out, dry = TRUE) {
  Sys.setenv(PATH=paste0("/usr/local/opt/libxml2/bin:/Users/shansabri/miniconda3/bin:/Library/Frameworks/Python.framework/Versions/3.4/bin:",
                         "/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/ncbi/blast/bin:/Library/TeX/texbin:/usr/local/MacGPG2/bin:",
                         "/opt/X11/bin:/Users/shansabri/bin:/Users/shansabri/edirect:/Users/shansabri/homer/bin/:/Users/shansabri/bin:/Users/shansabri/edirect"))
  
  cmd <- sprintf("findMotifsGenome.pl %s %s %s -p 7 -size given", file, fa, out) # -bg %s, bg
  
  if(dry == TRUE) {
    return(cmd)
  } else {
    try(system(cmd))  
    return(TRUE)
  }
}


CMDS <- lapply(seq_along(programs), function(x){
  prog <- programs[[x]]
  prog_id <- names(programs)[x]; message(prog_id)
  prog_info <- get_gene_info(prog) 
  prog_info <- reduce_to_tss(prog_info)
  prog_info <- extend(prog_info, upstream = 500, downstream = 500)
  
  outdir <- file.path(getwd(), prog_id); dir.create(outdir, showWarnings = FALSE)
  
  
  if(grepl("1111", prog_id)) {
    states <- read.table("~/Dropbox/ErnstLab/deconvolution_old/CHIP/chromHMM/NODOX_18_dense.bed", header = TRUE, row.names = NULL)
    message("Using states from NODOX_18_dense.bed.....")
  } else if(grepl("3333", prog_id)) {
    states <- read.table("~/Dropbox/ErnstLab/deconvolution_old/CHIP/chromHMM/ES_18_dense.bed", header = TRUE, row.names = NULL)
    message("Using states from ES_18_dense.bed.....")
  } else {
    message("PROG_ID does not contain 1111 or 3333")
  }
  
  states <- subset(states, description.. == "1_PromA")
  states <- GenomicRanges::GRanges(
    seqnames = Rle(states[[1]]),
    ranges = IRanges(start = states[[2]], end = states[[3]], names = seq(1, nrow(states))),
    state = toupper(states[[4]])
  )

  # intersect states with program info
  hits <- findOverlaps(query = states, subject = prog_info, type = "any", select = "all")
  state=states$state[hits@from]
  state_start=start(states)[hits@from]
  state_end=end(states)[hits@from]
  gene_symbol=prog_info$symbol[hits@to]
  gene_chr=seqnames(prog_info)[hits@to]
  gene_start=ifelse(strand(prog_info)[hits@to] == '+', start(prog_info)[hits@to], end(prog_info)[hits@to])
  gene_end=ifelse(strand(prog_info)[hits@to] == '+', end(prog_info)[hits@to], start(prog_info)[hits@to])
  gene_strand=strand(prog_info)[hits@to]
  overlaps <- data.frame(gene_symbol=gene_symbol,
                         gene_chr=gene_chr,
                         gene_start=gene_start,
                         gene_end=gene_end,
                         gene_strand=gene_strand,
                         state=state,
                         state_start =  state_start,
                         state_end = state_end)
  
  homer_input <- file.path(outdir, "state1_tss_overlap.bed")
  write.table(overlaps[,c("gene_chr", "state_start", "state_end")], homer_input, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

  cmd <- run_homer(file = homer_input, 
                   fa = "~/Dropbox/Binhouse/mm9.fa.gz",
                   bg = NULL,
                   out = outdir, 
                   dry = TRUE)
  return(cmd)
})






# #### ONLY TSS
# cmds <- lapply(seq_along(programs), function(n){
#   # n <- 1
#   id <- names(programs)[n]
#   outdir <- file.path(getwd(), id)
#   dir.create(outdir, showWarnings = FALSE)
#   for_homer <- convert_genes(programs[[n]])
#   loss <- round((length(unique(as.vector(programs[[n]]$Genes))) - nrow(for_homer)) / length(unique(as.vector(programs[[n]]$Genes))), 3)
#   message(paste(id, "Biomart lookup loss percent:", loss))
#   # write.table(for_homer, "~/.Trash/tmp.txt", quote = FALSE, row.names = FALSE)
#   homer_input_tss <- file.path(outdir, "tss.txt")
#   write.table(for_homer, homer_input_tss, quote = FALSE, row.names = FALSE)
#   cmd <- run_homer(file = homer_input_tss, out = outdir, dry = TRUE)
#   return(cmd)
# })