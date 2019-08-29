### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Shawn_Tan/2019-08-19/analysis/find-diffexp-genes/")



### LIBS
pacman::p_load(ggplot2, ggrepel, DESeq2, dplyr, tidyr)
theme_set(theme_bw(base_size = 14) + 
            theme(panel.grid = element_blank(), 
                  legend.justification = "top"))


### DATA
count_dir <- '~/Dropbox/PlathLab/Analyses/Shawn_Tan/2019-08-19/counts/'
files <- list.files(count_dir, pattern="*.counts", full.names=TRUE)
data <- do.call(cbind.data.frame, pbapply::pblapply(files, FUN=function(files){read.table(files,header=F, row.names=1, sep="\t")}))
data <- head(data, -5) 
colnames(data) <- gsub("\\.counts", "", basename(files))
colnames(data) <- gsub("-", "_", colnames(data))
colnames(data) <- sapply(strsplit(colnames(data), "_S", fixed=TRUE), function(x) (x[1]))
mapping <- readRDS("~/Dropbox/PlathLab/Analyses/Shawn_Tan/2019-08-19/analysis/data/mapping.rds")


### COMP1 
# samples <- c("1_day_DOX", "1_day_no_DOX")

### COMP2
samples <- c("2_days_DOX", "2_days_no_DOX", "1_day_DOX_1_day_no_DOX")

### COMP3
# samples <- c("3_days_DOX", "3_days_no_DOX", "1_day_DOX_2_days_no_DOX", "2_days_DOX_1_day_no_DOX")

### COMP4
# samples <- c("4_days_DOX", "4_days_no_DOX", "2_days_DOX_2_days_no_DOX", "3_days_DOX_1_day_no_DOX")

### COMP5
# samples <- c("5_days_DOX", "5_days_no_DOX", "3_days_DOX_2_days_no_DOX")

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
rld <- rld[rowMeans(rld) >= quantile(rowMeans(rld), 0.75), ]




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





















outdir <- file.path(getwd(), paste(samp1, samp2, sep = "-vs-"))
dir.create(outdir, showWarnings = FALSE)

data <- counts[, c("chromosome_name", samp1, samp2)]
names(data) <- c("chromosome_name", "A", "B")
# data <- subset(data, chromosome_name == "X")
data$FC <- data$A - data$B
data$Gene <- row.names(data)
data$chr_type <- ifelse(data$chromosome_name %in% c("X"), "chrX", "chrAUTOSOMAL")


## FC DISTRIBUTION 
# ggplot(data, aes(FC)) +
#   geom_density(aes(y = ..scaled..), fill = "deeppink4", alpha = 0.5) + 
#   labs(x = "Rlog2 Fold Change", y = "Scaled Density", 
#        title = paste(samp1, samp2, sep = " vs. "), 
#        subtitle = paste(nrow(data), "X-linked Genes")) +
#   ggsave(file.path(outdir, "FOLD-CHANGE-DISTRIBUTION.png"), height = 5, width = 9)

chr_type_tbl <- table(data$chr_type)
ggplot(data, aes(FC)) +
  facet_wrap(~ chr_type, ncol = 1, scales = "free") +
  geom_density(aes(y = ..scaled..), fill = "deeppink4", alpha = 0.5) + 
  labs(x = "Rlog2 Fold Change", y = "Scaled Density", 
       title = paste(samp1, samp2, sep = " vs. "), 
       subtitle = sprintf("%s X-linked Genes\n%s Autosomal Genes", chr_type_tbl[2], chr_type_tbl[1])) +
  ggsave(file.path(outdir, "FOLD-CHANGE-DISTRIBUTION.png"), height = 7, width = 9)



## DEFINED THRESHOLDS
fc_thresold <- 1
cols <- c("Significant" = "red", "Not Significant" = "grey")

data$Significant <- ifelse(abs(data$FC) >= 1, "Significant", "Not Significant")
table(data$Significant, data$chr_type)

# ggplot(data, aes(FC)) +
#   geom_density(aes(y = ..scaled..), fill = "deeppink4", alpha = 0.5) + 
#   geom_vline(xintercept = c(-1*fc_thresold, fc_thresold), linetype="dashed", color = "deeppink4", size=1) +
#   labs(x = "Rlog2 Fold Change", y = "Scaled Density", 
#        title = paste(samp1, samp2, sep = " vs. "), 
#        subtitle = paste(nrow(data), "X-linked Genes; abs(Rlog2 FC) critera = ", fc_thresold)) +
#   ggsave(file.path(outdir, "FOLD-CHANGE-DISTRIBUTION-WITH-FC-CUTOFF.png"), height = 5, width = 9)

ggplot(data, aes(FC)) +
  facet_wrap(~ chr_type, ncol = 1, scales = "free") +
  geom_density(aes(y = ..scaled..), fill = "deeppink4", alpha = 0.5) + 
  geom_vline(xintercept = c(-1*fc_thresold, fc_thresold), linetype="dashed", color = "deeppink4", size=1) +
  labs(x = "Rlog2 Fold Change", y = "Scaled Density", 
       title = paste(samp1, samp2, sep = " vs. "), 
       subtitle = sprintf("%s X-linked Genes\n%s Autosomal Genes\nabs(Rlog2 FC) critera of %s", 
                          chr_type_tbl[2], chr_type_tbl[1], fc_thresold)) +
  ggsave(file.path(outdir, "FOLD-CHANGE-DISTRIBUTION-WITH-FC-CUTOFF.png"), height = 7, width = 9)


ggplot(data, aes(A, B)) +
  facet_wrap(~ chr_type, nrow = 1, scales = "free") +
  geom_point(aes(colour = Significant)) +
  # geom_text_repel(data = subset(data, Significant == "Significant"), aes(label = Gene), 
  geom_text_repel(data = subset(data, Significant == "Significant" & chr_type == "chrX"), aes(label = Gene), 
                  size = 2, segment.alpha = 0.4, segment.size = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed") +
  scale_colour_manual(values = cols) +
  labs(x = paste0("Rlog2(", samp1, ")"), y = paste0("Rlog2(", samp2, ")"), 
       title = paste(samp1, samp2, sep = " vs. "), 
       subtitle = paste0("abs(Rlog2 FC) critera = ", fc_thresold, "; ", nrow(subset(data, Significant == "Significant")), " differential genes")) +
  ggsave(file.path(outdir, "SCATTER-DIFF-EXP.png"), height = 7, width = 15)




### EXPORT 
names(data)[names(data) == "A"] = samp1
names(data)[names(data) == "B"] = samp2
saveRDS(data, file = file.path(outdir, sprintf("%s-vs-%s.rds", samp1, samp2)), compress = TRUE)
