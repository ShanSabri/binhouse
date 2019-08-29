### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/dir/")



### LIBS
pacman::p_load(ggplot2, reshape2, DESeq2, biomaRt, ggsignif, dplyr, tidyr)
theme_set(theme_bw(base_size = 14) + 
            theme(panel.grid = element_blank(), 
                  legend.justification = "top"))



### FXNS
gatherpairs <- function(data, ..., 
                        xkey = '.xkey', xvalue = '.xvalue',
                        ykey = '.ykey', yvalue = '.yvalue',
                        na.rm = FALSE, convert = FALSE, factor_key = FALSE) {
  vars <- quos(...)
  xkey <- enquo(xkey)
  xvalue <- enquo(xvalue)
  ykey <- enquo(ykey)
  yvalue <- enquo(yvalue)
  
  data %>% {
    cbind(gather(., key = !!xkey, value = !!xvalue, !!!vars,
                 na.rm = na.rm, convert = convert, factor_key = factor_key),
          select(., !!!vars)) 
  } %>% gather(., key = !!ykey, value = !!yvalue, !!!vars,
               na.rm = na.rm, convert = convert, factor_key = factor_key)
}



### DATA
count_dir <- '~/Dropbox/PlathLab/Analyses/Shawn_Tan/2019-08-19/counts/'
files <- list.files(count_dir, pattern="*.counts", full.names=TRUE)
data <- do.call(cbind.data.frame, pbapply::pblapply(files, FUN=function(files){read.table(files,header=F, row.names=1, sep="\t")}))
data <- head(data, -5) 
colnames(data) <- gsub("\\.counts", "", basename(files))
colnames(data) <- gsub("-", "_", colnames(data))
colnames(data) <- sapply(strsplit(colnames(data), "_S", fixed=TRUE), function(x) (x[1]))
# write.table(data, file="counts.tsv", sep="\t", row.names=T, col.names=NA)  

# mart <- useDataset(
#   "mmusculus_gene_ensembl",
#   useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "may2012.archive.ensembl.org")
# )
# mapping <- getBM(
#   filters = "mgi_symbol",
#   attributes = c("mgi_symbol", "chromosome_name", "start_position", "end_position", "strand"),
#   values = row.names(data),
#   mart = mart
# )
# saveRDS(mapping, compress = TRUE, "~/Dropbox/PlathLab/Analyses/Shawn_Tan/2019-08-19/analysis/data/mapping.rds")
mapping <- readRDS("~/Dropbox/PlathLab/Analyses/Shawn_Tan/2019-08-19/analysis/data/mapping.rds")


### COMP1 
# samples <- c("1_day_DOX", "1_day_no_DOX")

### COMP2
# samples <- c("2_days_DOX", "2_days_no_DOX", "1_day_DOX_1_day_no_DOX")

### COMP3
# samples <- c("3_days_DOX", "3_days_no_DOX", "1_day_DOX_2_days_no_DOX", "2_days_DOX_1_day_no_DOX")

### COMP4
# samples <- c("4_days_DOX", "4_days_no_DOX", "2_days_DOX_2_days_no_DOX", "3_days_DOX_1_day_no_DOX")

### COMP5
samples <- c("5_days_DOX", "5_days_no_DOX", "3_days_DOX_2_days_no_DOX")

counts <- data[, samples]
metadata <- data.frame(row.names = names(counts),
                       ID = names(counts), 
                       CONDITION = sapply(strsplit(names(counts), "_", fixed=TRUE), function(x) (paste(x[1], x[2], sep = "_"))))

### FILTER 
idx <- which(rowSums(counts) == 0)
counts <- counts[-idx,]
keep <- rowSums(counts>3) >= 2 ; table(keep)
counts <- counts[keep, , drop=F]




### DESEQ NORMALIZATION 
# dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design= ~ CONDITION)
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design= ~ 1)
dds <- DESeq(dds)
rld <- as.data.frame(assay(rlog(dds, blind=FALSE)))




### FILTER NORMALIZED COUNTS FOR HIGH EXPRESSED GENES
# rld <- rld[rowMeans(rld) >= quantile(rowMeans(rld), 0.75), ]




### TAG GENE WITH CHR
rld <- merge(mapping, rld, by.x = "mgi_symbol", by.y = 0)
rld <- subset(rld, chromosome_name %in% c(seq(1, 19), "X"))
# row.names(rld) <- rld$mgi_symbol; rld$mgi_symbol <- NULL
row.names(rld) <- make.names(rld$mgi_symbol, unique = TRUE)
# write.table(rld, file=gzfile("~/Dropbox/PlathLab/Analyses/Shawn_Tan/2019-07-11/data/COUNTS-RLD-NORM.tsv.gz"),
#             sep="\t", col.names=NA, row.names=TRUE, quote = FALSE)
# saveRDS(rld, compress = TRUE, "~/Dropbox/PlathLab/Analyses/Shawn_Tan/2019-07-11/data/COUNTS-RLD-NORM.rds")



### PLOT
data <- rld
data$chr_id <- ifelse(data$chromosome_name %in% seq(1, 19), "Autosomal", "X")
data_m <- melt(data, id.vars = c("chromosome_name", "chr_id", "start_position", "end_position", "strand", "mgi_symbol"))
data_m$DOX <- ifelse(grepl("_no_DOX", data_m$variable), "No DOX", "DOX")
n_genes <- table(data$chr_id)



### SCATTER PLOTS
data %>% 
  gatherpairs(samples) %>% 
  group_by(chr_id, .xkey, .ykey) %>% 
  summarise(r = paste(" P_r =", round(cor(.xvalue, .yvalue), 5))) -> cors

A <- data %>% 
  gatherpairs(samples) %>% 
  filter(chr_id == "X") %>% {
    ggplot(., aes(x = .xvalue, y = .yvalue)) +
      facet_wrap(.xkey ~ .ykey, ncol = length(unique(.$.ykey)), scales = 'free', labeller = label_both) +
      geom_point(alpha = 0.5, size = 0.5) + 
      geom_text(data = subset(cors, chr_id == "X"), aes(label = r, x = -Inf, y = Inf), hjust = 0, vjust = 1.5) +
      geom_smooth(method = 'lm') +
      scale_color_brewer(type = 'qual') +
      labs(title = "chrX") 
      # ggsave("PAIRWISE-SCATTER.png", height = 10, width = 10)
    # ggsave("PAIRWISE-SCATTER-X.png", height = 15, width = 15)
  }

B <- data %>% 
  gatherpairs(samples) %>% 
  filter(chr_id == "Autosomal") %>% {
    ggplot(., aes(x = .xvalue, y = .yvalue)) +
      facet_wrap(.xkey ~ .ykey, ncol = length(unique(.$.ykey)), scales = 'free', labeller = label_both) +
      geom_point(alpha = 0.5, size = 0.5) + 
      geom_text(data = subset(cors, chr_id == "Autosomal"), aes(label = r, x = -Inf, y = Inf), hjust = 0, vjust = 1.5) +
      geom_smooth(method = 'lm') +
      scale_color_brewer(type = 'qual') +
      labs(title = "chrAutosomal") 
      # ggsave("PAIRWISE-SCATTER.png", height = 10, width = 10)
      # ggsave("PAIRWISE-SCATTER-AUTO.png", height = 15, width = 15)
  }

PLOT <- cowplot::plot_grid(A, B) +
  ggsave("PAIRWISE-SCATTER.png", height = 15, width = 30)



# ### COMP1
# ggplot(data_m, aes(x = variable, y = value, group = variable, fill = DOX)) +
#   facet_wrap(~ chr_id, nrow = 1) +
#   geom_boxplot(alpha = 0.2, position = "dodge2") +
#   geom_signif(comparisons = list(samples), map_signif_level=FALSE) +
#   labs(x = "", y = "rlog(counts)",
#        title = "Distribution of chrAuto and chrX normalized gene expression",
#        subtitle = paste0(n_genes[1], " Autosomal genes & ", n_genes[2], " X-linked genes")) +
#   # ylim(0, 21) +
#   ylim(10, 21) +
#   theme_bw(base_size = 16) +
#   theme(strip.text = element_text(size = 16),
#         axis.text.x = element_text(angle = 60, hjust = 1),
#         panel.grid = element_blank(),
#         legend.position = "none") +
#   ggsave("BOX-DIST-AUTO-VS-X-SIGNIF-PVAL.png", height = 8, width = 10)


# head(data_m)
# head(data)
# ggplot(data, aes(x = `1_day_DOX`, y = `1_day_no_DOX`)) +
#   facet_wrap(~ chr_id, nrow = 1) +
#   geom_point(alpha = 0.2)
# 
# # wilcox.test(subset(data, chr_id == "Autosomal")$`1_day_DOX`, subset(data, chr_id == "Autosomal")$`1_day_no_DOX`)$p.val
# # wilcox.test(subset(data, chr_id == "X")$`1_day_DOX`, subset(data, chr_id == "X")$`1_day_no_DOX`)$p.val


# ### COMP2
# ggplot(data_m, aes(x = variable, y = value, group = variable, fill = DOX)) +
#   facet_wrap(~ chr_id, nrow = 1) +
#   geom_boxplot(alpha = 0.2, position = "dodge2") +
#   geom_signif(comparisons = combn(samples, 2, list), map_signif_level=FALSE, y_position = c(19,20,19)) +
#   labs(x = "", y = "rlog(counts)",
#        title = "Distribution of chrAuto and chrX normalized gene expression",
#        subtitle = paste0(n_genes[1], " Autosomal genes & ", n_genes[2], " X-linked genes")) +
#   # ylim(0, 21) +
#   ylim(10, 21) +
#   theme_bw(base_size = 16) +
#   theme(strip.text = element_text(size = 16),
#         axis.text.x = element_text(angle = 60, hjust = 1),
#         panel.grid = element_blank(),
#         legend.position = "none") +
#   ggsave("BOX-DIST-AUTO-VS-X-SIGNIF-PVAL.png", height = 8, width = 10)

# ### COMP3
# ggplot(data_m, aes(x = variable, y = value, group = variable, fill = DOX)) +
#   facet_wrap(~ chr_id, nrow = 1) +
#   geom_boxplot(alpha = 0.2, position = "dodge2") +
#   geom_signif(comparisons = combn(samples, 2, list), map_signif_level=FALSE, y_position = c(19,20,21,22,23,24)+1) +
#   labs(x = "", y = "rlog(counts)",
#        title = "Distribution of chrAuto and chrX normalized gene expression",
#        subtitle = paste0(n_genes[1], " Autosomal genes & ", n_genes[2], " X-linked genes")) +
#   # ylim(0, 25) +
#   ylim(10, 25) +
#   theme_bw(base_size = 16) +
#   theme(strip.text = element_text(size = 16),
#         axis.text.x = element_text(angle = 60, hjust = 1),
#         panel.grid = element_blank(),
#         legend.position = "none") +
#   ggsave("BOX-DIST-AUTO-VS-X-SIGNIF-PVAL.png", height = 8, width = 10)

# ### COMP4
# ggplot(data_m, aes(x = variable, y = value, group = variable, fill = DOX)) +
#   facet_wrap(~ chr_id, nrow = 1) +
#   geom_boxplot(alpha = 0.2, position = "dodge2") +
#   geom_signif(comparisons = combn(samples, 2, list), map_signif_level=FALSE, y_position = c(19,20,21,22,23,24)+1) +
#   labs(x = "", y = "rlog(counts)",
#        title = "Distribution of chrAuto and chrX normalized gene expression",
#        subtitle = paste0(n_genes[1], " Autosomal genes & ", n_genes[2], " X-linked genes")) +
#   # ylim(0, 25) +
#   ylim(10, 25) +
#   theme_bw(base_size = 16) +
#   theme(strip.text = element_text(size = 16),
#         axis.text.x = element_text(angle = 60, hjust = 1),
#         panel.grid = element_blank(),
#         legend.position = "none") +
#   ggsave("BOX-DIST-AUTO-VS-X-SIGNIF-PVAL.png", height = 8, width = 10)

# ### COMP5
# ggplot(data_m, aes(x = variable, y = value, group = variable, fill = DOX)) +
#   facet_wrap(~ chr_id, nrow = 1) +
#   geom_boxplot(alpha = 0.2, position = "dodge2") +
#   geom_signif(comparisons = combn(samples, 2, list), map_signif_level=FALSE, y_position = c(19,20,19)+1) +
#   labs(x = "", y = "rlog(counts)",
#        title = "Distribution of chrAuto and chrX normalized gene expression",
#        subtitle = paste0(n_genes[1], " Autosomal genes & ", n_genes[2], " X-linked genes")) +
#   # ylim(0, 22) +
#   ylim(10, 22) +
#   theme_bw(base_size = 16) +
#   theme(strip.text = element_text(size = 16),
#         axis.text.x = element_text(angle = 60, hjust = 1),
#         panel.grid = element_blank(),
#         legend.position = "none") +
#   ggsave("BOX-DIST-AUTO-VS-X-SIGNIF-PVAL.png", height = 8, width = 10)
