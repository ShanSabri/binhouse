### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/Drop12/")



### LIBS
pacman::p_load(ggplot2, reshape2, RColorBrewer)
theme_set(theme_bw(base_size = 14) + 
            theme(panel.grid = element_blank(), 
                  legend.justification = "top"))



### DATA 
data_dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data"
raw <- readRDS(paste(data_dir, "DROP12-RAW-FILTERED-DGE.rds", sep = "/"))
norm <- log((sweep(raw, 2, apply(raw, 2, function(x) sum(x)), "/")*10000)+1)
pheno <- readRDS(paste(data_dir, "DROP12-PHENODATA.rds", sep = "/"))
identical(names(norm), row.names(pheno))



### X-LINKED GENES VS AUTO
generate_gene_annotation_biomart <- function(g) {
  mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl", host = "may2012.archive.ensembl.org"))
  mapping <- getBM(filters = "mgi_symbol", attributes = c("mgi_symbol", "chromosome_name"), values = g, mart = mart)
  mapping <- subset(mapping, chromosome_name %in% c(seq(1, 19), "X", "Y"))
  row.names(mapping) <- toupper(mapping$mgi_symbol);  mapping$mgi_symbol <- NULL
  return(mapping)
}
range_norm <- function(x){(x-min(x))/(max(x)-min(x))}
gene_chr_info <- generate_gene_annotation_biomart(row.names(norm))
pheno$chr_auto <- range_norm(colMeans(norm[row.names(subset(gene_chr_info, chromosome_name %in% seq(1, 19))),]))
pheno$chr_x <- range_norm(colMeans(norm[row.names(subset(gene_chr_info, chromosome_name %in% c("X"))),]))
pheno$x_minus_auto <- pheno$chr_x - pheno$chr_auto

ggplot(pheno[with(pheno, order(chr_auto)),], aes(tSNE_1, tSNE_2)) +
  geom_point(data=pheno[,c('tSNE_1', 'tSNE_2')], colour="grey", size=0.5) +  
  geom_point(aes(colour=chr_auto), size=0.5) +
  scale_color_distiller(palette = "RdPu", direction = 1) +
  theme(strip.text.x = element_text(size = 18)) +
  theme(axis.title=element_text(size=14), axis.text = element_text(size=10)) +
  ggsave(paste(getwd(), "TSNE-AUTO-CHR-FACET-RANGE-NORM.png", sep="/"), height=6, width=7)

ggplot(pheno[with(pheno, order(chr_x)),], aes(tSNE_1, tSNE_2)) +
  geom_point(data=pheno[,c('tSNE_1', 'tSNE_2')], colour="grey", size=0.5) +  
  geom_point(aes(colour=chr_x), size=0.5) +
  scale_color_distiller(palette = "RdPu", direction = 1) +
  theme(strip.text.x = element_text(size = 18)) +
  theme(axis.title=element_text(size=14), axis.text = element_text(size=10)) +
  ggsave(paste(getwd(), "TSNE-X-CHR-FACET-RANGE-NORM.png", sep="/"), height=6, width=7)

ggplot(pheno[with(pheno, order(x_minus_auto)),], aes(tSNE_1, tSNE_2)) +
  geom_point(data=pheno[,c('tSNE_1', 'tSNE_2')], colour="grey", size=0.5) +  
  geom_point(aes(colour=x_minus_auto), size=0.5) +
  # scale_color_distiller(palette = "PiYG", direction = 1) +
  scale_colour_gradient2() +
  theme(strip.text.x = element_text(size = 18)) +
  theme(axis.title=element_text(size=14), axis.text = element_text(size=10)) +
  ggsave(paste(getwd(), "TSNE-X-AUTO-CHR-FACET.png", sep="/"), height=6, width=7)



### XIST
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  substr(x, 2, nchar(x)) <- tolower(substr(x, 2, nchar(x)))
  return(x)
}

plot_single_gene <- function(g) {
  message(g)
  gene <- cbind.data.frame(pheno, t(norm[g, row.names(pheno)]))
  gene <- gene[with(gene, order(gene[[g]])), ]
  p <- ggplot(gene, aes(tSNE_1, tSNE_2)) +
    geom_point(aes_string(colour = g)) +
    scale_colour_gradientn(colours = c("grey55", "grey85", "red", "purple"), na.value = "grey80", "") +
    labs(title = firstup(g)) +
    theme(
      axis.title = element_text(size = 30),
      axis.text = element_text(size = 20),
      plot.title = element_text(hjust = 0.5, size = 38, face = "italic"),
      legend.text = element_text(size = 24),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.width = unit(2, "cm")
    ) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    ggsave(paste(getwd(), paste0("TSNE-", g, ".png"), sep = "/"), height = 8, width = 7)
  return(p)
}
plot_single_gene("XIST")
