################################################
## OPTS
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/factor-binding-analysis-kc/"); getwd()




################################################
## LIBS
pacman::p_load(
  biomaRt, GenomicRanges,
  pheatmap, RColorBrewer, 
  tidyverse, reshape2, rwantshue
)




################################################
## FXNS
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

extend <- function(x, upstream=0, downstream=0) {
  if (any(strand(x) == "*")){ warning("'*' ranges were treated as '+'")}
  on_plus <- strand(x) == "+" | strand(x) == "*"
  new_start <- start(x) - ifelse(on_plus, upstream, downstream)
  new_end <- end(x) + ifelse(on_plus, downstream, upstream)
  ranges(x) <- IRanges(new_start, new_end)
  # trim(x)
  return(x)
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

# zero_one <- function(x){
#   return( (x-min(x))/(max(x)-min(x)))
# }
# 
# min_max <- function(data, min, max) {
#   data2 <- data
#   data2[data2 > max] <- max
#   data2[data2 < min] <- min
#   return(data2)
# }


################################################
## BINDING DATA FOR ALL FACTORS AND PATTERNS
dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/kostas_factor_binding_patterns"
beds <- list.files(dir, pattern = "*.bed", full.names = TRUE)
binding <- lapply(beds, bed_to_granges)
names(binding) <- basename(beds)
saveRDS(binding, file.path(dir, "BINDING-PATTERNS.rds"), compress = TRUE)




################################################
## GET PROGRAM GENE SET INFO
programs <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/Drop10.Aggregated.Programs.Gene.Master.List.rds")
programs_feat <- generate_gene_annotation_biomart(as.vector(unique(programs$Genes)))
programs_feat <- reduce_to_tss(programs_feat)
programs <- merge(programs, as.data.frame(programs_feat), by.x = "Genes", by.y = "symbol", all.y = TRUE)
programs <- split(programs, f = programs$Set)
programs <- lapply(programs, function(x) with(x, GRanges(seqnames, IRanges(start, end), gene_id = x$Genes)))
annot <- data.frame(size = sapply(programs, length))




################################################
## INTERSECT - BACKGROUND
all_genes <- read.table("~/Dropbox/PlathLab/Code/Jubri/v2/all-genes.txt")
names(all_genes)[names(all_genes) == "V1"] <- "Genes"
all_genes <- generate_gene_annotation_biomart(g = as.vector(all_genes$Genes))
all_genes <- reduce_to_tss(all_genes)
all_genes <- extend(all_genes, upstream = 20000, downstream = 20000)
total_num_genes <- length(all_genes)
genes_bound_for_each_pattern_bg <- data.frame(tmp = sapply(seq_along(binding), function(y){
  binding_pattern <- sub("(.*?_.*?)_.*", "\\1", names(binding)[y])
  hits <- findOverlaps(query = binding[[y]], subject = all_genes, type = "any", select = "all")
  num_genes_bound <- length(unique(as.vector(all_genes$symbol[hits@to])))
  names(num_genes_bound) <- binding_pattern
  return(num_genes_bound)
}))
names(genes_bound_for_each_pattern_bg) <- "bg"
genes_bound_for_each_pattern_bg <- genes_bound_for_each_pattern_bg / total_num_genes


# PLOT
pheatmap(t(genes_bound_for_each_pattern_bg), 
         display_numbers = TRUE,
         cluster_rows = FALSE, 
         number_format = "%.3f",
         main = "BACKGROUND",
         color = colorRampPalette(brewer.pal(n = 7, name = "Blues"))(100),
         filename = paste(getwd(), "HEATMAP-BACKGROUND-RAW.png", sep="/"),
         height=3, width=15)




################################################
## INTERSECT - FOREGROUND
results <- lapply(seq_along(programs), function(x){
  id <- names(programs)[x]
  message(paste(x, id))
  
  x_ext <- extend(programs[[x]], upstream = 20000, downstream = 20000)
  
  genes_bound_for_each_pattern <- data.frame(tmp = sapply(seq_along(binding), function(y){
    binding_pattern <- sub("(.*?_.*?)_.*", "\\1", names(binding)[y])
    hits <- findOverlaps(query = binding[[y]], subject = x_ext, type = "any", select = "all")
    num_genes_bound <- length(unique(as.vector(x_ext$gene_id[hits@to])))
    names(num_genes_bound) <- binding_pattern
    return(num_genes_bound)
  }))
  names(genes_bound_for_each_pattern) <- id 
  num_genes_in_prog <- length(programs[[x]])
  genes_bound_for_each_pattern <- genes_bound_for_each_pattern / num_genes_in_prog
  return(genes_bound_for_each_pattern)
})
results_df <- do.call(cbind.data.frame, results)

names(results) <-  names(programs)
saveRDS(results, paste(getwd(), "PROGRAM-INTERSECT-RESULTS.rds", sep="/"), compress = TRUE)

# PLOT
pheatmap(results_df, 
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "FOREGROUND",
         annotation_col = annot, 
         annotation_names_col = FALSE,
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         filename = paste(getwd(), "HEATMAP-FOREGROUND-RAW.png", sep="/"),
         height=8, width=15)

pheatmap(results_df, 
         display_numbers = TRUE,
         number_format = "%.2f",
         cluster_rows = FALSE,
         annotation_col = annot, 
         annotation_names_col = FALSE,
         main = "FOREGROUND",
         color = colorRampPalette(brewer.pal(n = 7, name = "Reds"))(100),
         filename = paste(getwd(), "HEATMAP-FOREGROUND-RAW-ORDERED.png", sep="/"),
         height=8, width=15)




################################################
## COMPUTE ENRICHMENT
enrichment <- do.call(cbind.data.frame, lapply(results, function(x) x / genes_bound_for_each_pattern_bg))
write.table(enrichment, file = file.path(getwd(), "ENRICHMENT-PCT.txt"), sep = "\t", col.names = NA)

# PLOT
pheatmap(enrichment, 
         display_numbers = TRUE,
         annotation_col = annot, 
         annotation_names_col = FALSE,
         number_format = "%.2f",
         main = "ENRICHMENT",
         color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100),
         filename = paste(getwd(), "HEATMAP-ENRICHMENT-RAW.png", sep="/"),
         height=8, width=15)

programs <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/Drop10.Aggregated.Programs.Gene.Master.List.rds")
order_of_progs <- unique(programs[, c("Order", "Set")])
pheatmap(enrichment[, as.vector(order_of_progs$Set)], 
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         display_numbers = TRUE,
         annotation_col = annot, 
         annotation_names_col = FALSE,
         number_format = "%.2f",
         main = "ENRICHMENT",
         color = colorRampPalette(brewer.pal(n = 7, name = "Greens"))(100),
         filename = paste(getwd(), "HEATMAP-ENRICHMENT-RAW-ORDERED.png", sep="/"),
         height=8, width=15)