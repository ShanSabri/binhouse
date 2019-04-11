## OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Tom_Allison/2019-03-20/")
pacman::p_load(GenomicAlignments, msa)




## BAMS
dir <- file.path(getwd(), "bams")
bams <- list.files(dir, pattern = "*.bam", full.names = TRUE, recursive = TRUE)
bams <- bams[!grepl(".bai", bams)]
samples <- unique(sapply(strsplit(tools::file_path_sans_ext(basename(bams)), "_", fixed=TRUE), function(x) (x[1])))



## UNIQ
chrA <- "OE_PLASMID"
chrB <- "ASCL1_BACKBONE"
chrC <- "DXL2_BACKBONE"

sapply(samples, function(x){

  # x = "N701"
  
  chrA_read_alnmt <- names(readGAlignments(file.path(getwd(), "bams", paste0(x, "_", chrA, ".bam")), use.names=TRUE))
  chrB_read_alnmt <- names(readGAlignments(file.path(getwd(), "bams", paste0(x, "_", chrB, ".bam")), use.names=TRUE))
  chrC_read_alnmt <- names(readGAlignments(file.path(getwd(), "bams", paste0(x, "_", chrC, ".bam")), use.names=TRUE))
  
  chrA_read_alnmt[1:5]
  chrB_read_alnmt[1:5]
  chrC_read_alnmt[1:5]
  
  shared_chrA_w_chrsBC <- min(0, intersect(chrA_read_alnmt, c(chrB_read_alnmt, chrC_read_alnmt)))
  shared_chrB_w_chrsAC <- min(0, intersect(chrB_read_alnmt, c(chrA_read_alnmt, chrC_read_alnmt)))
  shared_chrC_w_chrsAB <- min(0, intersect(chrC_read_alnmt, c(chrA_read_alnmt, chrB_read_alnmt)))
  
  message(sprintf("%s\n\tShared between A and BC: %s\n\tShared between B and AC: %s\n\tShared between C and AB: %s", 
          x, shared_chrA_w_chrsBC, shared_chrB_w_chrsAC, shared_chrC_w_chrsAB))
})


## MSA TO FIND OVERLAPPING REGIONS
fq <- file.path(getwd(), "new_chrs.fa")
seq <- readDNAStringSet(fq)
alignment <- msa(seq, "ClustalOmega", verbose = TRUE)
print(alignment, showNames=FALSE, show="complete")
tex_file <- file.path(getwd(), "alignment.tex")
# msaPrettyPrint(alignment, output="tex", showNames="none", showLogo="none", askForOverwrite=FALSE, file = tex_file)
# tools::texi2pdf(tex_file, clean=TRUE)

msaPrettyPrint(alignment, output="pdf", showNames="none", showLogo="none", askForOverwrite=FALSE, file = gsub("tex", "pdf", tex_file))

chunkSize <- 2500
for (start in seq(1, ncol(alignment), by=chunkSize)) {
  end <- min(start + chunkSize - 1, ncol(alignment))
  alnPart <- DNAMultipleAlignment(subseq(unmasked(alignment), start, end))
  msaPrettyPrint(x=alnPart, output="pdf", subset=NULL, file=file.path(getwd(), "msa", paste0("aln_", start, "-", end, ".pdf")))
}




