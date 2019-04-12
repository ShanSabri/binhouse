## OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Tom_Allison/2019-03-20/")
pacman::p_load(GenomicAlignments, rtracklayer)



## FILTER BED FILES BY READS IN BAM
samples <- paste0("N", c(701, 704, 705))
CATCH <- sapply(samples, function(x){
  message(sprintf("%s: Processing %s", Sys.time(), x))
  bam_file <- file.path(getwd(), "bams", "cleaned", paste0(x, "_OE_PLASMID.bam"))
  bed_file <- file.path(getwd(), "beds", paste0(x, "_hg19_OE_PLASMID.bed") )
  out_file <- file.path(getwd(), "beds", "cleaned", paste0(x, "_hg19_OE_PLASMID.bed"))
  # length(aln)
  # length(bed)
  # length(intersect(mcols(aln)$qname, mcols(bed)$name))
  
  bed <- import(bed_file, format="bed")
  aln <- readGAlignments(bam_file, param=ScanBamParam(what="qname"))
  bed_filt <- bed[mcols(bed)$name %in% mcols(aln)$qname]
  export(bed_filt, format = "bed", con = out_file)
})
