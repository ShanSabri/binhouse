### OPT
rm(list = ls(all = TRUE))




### READ DATA
dge <- data.table::fread("~/Downloads/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  nThread = parallel::detectCores() - 1,
  showProgress = TRUE,
  verbose = TRUE
)




### REDUCED DATA TO GENE LEVEL
library(tidyverse)
## exploratory -- number of transcripts per gene
## note -- wtf why is there a gene with 100+ transcripts
# dge %>%
#   group_by(gene_id) %>%
#   summarise(num_transcripts = n()) %>%
#   ggplot(., aes(num_transcripts)) +
#     geom_density(aes(y = ..scaled..), fill = "blue", alpha = 0.5) +
#     scale_x_continuous(trans='log10') +
#     theme_minimal(base_size = 14)
dge %>%
  mutate(cleaned_gene_id = sapply(strsplit(gene_id, ".", fixed = T), function(x) x[1])) %>%
  select(-c(gene_id, transcript_id)) %>%
  group_by(cleaned_gene_id) %>%
  summarise_all(funs(mean)) -> dge_reduced
dge_reduced <- data.frame(dge_reduced)
row.names(dge_reduced) <- toupper(dge_reduced$cleaned_gene_id)
dge_reduced$cleaned_gene_id <- NULL
rm(dge)
gc() # free up memory




### TAG DATA WITH GENE ID
library(biomaRt)
mart <- useDataset(
  "hsapiens_gene_ensembl",
  useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl")
)
mapping <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = row.names(dge_reduced),
  mart = mart
)
dge <- merge(mapping, dge_reduced, by.x = "ensembl_gene_id", by.y = 0)
dge$ensembl_gene_id <- NULL
row.names(dge) <- make.names(toupper(dge$hgnc_symbol), unique = TRUE)
dge$hgnc_symbol <- NULL
dge[1:10, 1:9]




### EXPORT
saveRDS(dge, compress = TRUE, file = "~/Downloads/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.rds")




### IMPORT
rm(list = ls(all = TRUE))
dge <- readRDS("~/Downloads/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.rds")
dim(dge)
dge[1:10, 1:9]
