### OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Haibin_Xi/2019-07-24/OE_CHRS/")



### LIBS
pacman::p_load(GenomicRanges)



### FUNCTIONS
bed_to_granges <- function(file) {
  message(file)
  df <- read.table(gzfile(file), header = F, stringsAsFactors = F)
  header <- c("chr", "start", "end", "id", "score", "strand", 
              "read2", "read1", "cell_barcode", "umi")
  names(df) <- header[1:length(names(df))]
  if ("strand" %in% colnames(df)) {
    df$strand <- gsub(pattern = "[^+-]+", replacement = "*", x = df$strand)
  }
  gr <- with(df, GRanges(chr, IRanges(start, end), id = id, score = score, strand = strand, 
                         read2 = read2, read1 = read1, cell_barcode = cell_barcode, umi = umi))

  return(gr)
}

strsplits <- function(x, splits, ...){
  for (split in splits) {
    x <- unlist(strsplit(x, split, ...))
  }
  return(x[!x == ""]) 
}

### FILES
dir <- file.path(getwd(), "intermediates")
beds <- list.files(dir, pattern = "*.txt.gz", full.names = TRUE, recursive = TRUE)
to_filter <- lapply(beds, bed_to_granges)
names(to_filter) <- basename(beds)




### FILTERING CRITERA 
critera <- with(data.frame(chr = c("OE_COMMON", "OE_TF", "OE_RTTA"), 
                           start = c(1675, 1086, 379), 
                           end = c(2271, 1796, 1142)), 
     GRanges(chr, IRanges(start, end)))



### FILTER
filtered <- lapply(seq_along(to_filter), function(x){
  id <- names(to_filter)[x]; message(id)
  chr <- paste(strsplits(id, c("_", "."), fixed=TRUE)[[2]], strsplits(id, c("_", "."), fixed=TRUE)[[3]], sep = "_")
  gr <- to_filter[[x]]
  coords_to_filter <- critera[seqnames(critera) == chr, ]
  
  
  hits <- findOverlaps(query = gr, subject = coords_to_filter, type = "any", select = "all",  ignore.strand = TRUE)
  start <- start(gr)[hits@from]
  end <- end(gr)[hits@from]
  strand <- strand(gr)[hits@from]
  ids <- gr$id[hits@from]
  score <- gr$score[hits@from]
  read2 <- gr$read2[hits@from]
  read1 <- gr$read1[hits@from]
  cell_barcode <- gr$cell_barcode[hits@from]
  umi <- gr$umi[hits@from]

  to_ret <- data.frame(chr = chr, start = start, end = end, strand = strand, 
                       id = id, score = score, read2 = read2, read1 = read1, 
                       cell_barcode = cell_barcode, umi = umi)
    
  write.table(to_ret, file=gzfile(file.path(dir, "filtered", gsub(".txt.gz", "_FILTERED.txt.gz", id))), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  write.table(to_ret[,c(1:3)], file=file.path(dir, "filtered", gsub(".txt.gz", "_FILTERED.bed", id)), 
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
  return(to_ret)
})














