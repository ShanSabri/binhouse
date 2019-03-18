### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/paper_figs/figure_1/")




### LIBS
pacman::p_load(ggplot2, reshape2, biomaRt)
theme_set(theme_bw(base_size = 14) + 
            theme(panel.grid = element_blank(), 
                  legend.justification = "top"))
pal <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")




### DATA
data_dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data"
norm <- readRDS(paste(data_dir, "TC1-NORM-FILTERED-DGE.rds", sep = "/"))
pheno <- readRDS(paste(data_dir, "TC1-PHENODATA.rds", sep = "/"))
identical(names(norm), row.names(pheno))




### GO
# id <- "GO:0030509" # BMP signaling pathway
# id <- "GO:0007179" # TGF beta receptor signaling pathway
# id <- "GO:0007219" # Notch signaling pathway
# id <- "GO:0016055" # Wnt signaling pathway
# id <- "GO:0070371" # ERK1 and ERK2 cascade
id <- "GO:0004707" # MAP kinase activity
mapping <- getBM(filters = "go", 
                 attributes = c("mgi_symbol"), 
                 values = id, 
                 mart = useDataset(dataset="mmusculus_gene_ensembl", useMart("ensembl"))) 
mapping$mgi_symbol <- toupper(mapping$mgi_symbol)
paste(mapping$mgi_symbol, collapse = ", "); nrow(mapping)


### PLOT
range_norm <- function(x){
  (x-min(x)) / (max(x)-min(x))
}

plot_go <- function(pheno, norm, genes, title = NULL){
  message(title)
  
  md <- t(subset(norm[, row.names(pheno)], row.names(norm) %in% toupper(genes)))
  md <- cbind.data.frame(pheno, n_markers_exp = range_norm(rowMeans(md)))  

    
  ggplot(md[with(md, order(n_markers_exp)), ], aes(tSNE_1, tSNE_2)) +
  # ggplot(md[with(md, order(n_markers_exp)), ], aes(UMAP1, UMAP2)) +
    geom_point(aes(colour = n_markers_exp)) +
    scale_colour_gradientn(colours = c("grey70","grey50","lightgoldenrod","orange","red"), name = "") +
    labs(subtitle = "Normalized Avg Expression", title = title) +
    theme(
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, size = 32),
      plot.subtitle = element_text(hjust = 0.5, size = 20),
      legend.text = element_text(size = 24),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.width = unit(2, "cm")
    ) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)) 
    ggsave(paste(getwd(), paste0("TSNE-GO-",gsub(" ", "-", title),".png"), sep = "/"), height = 8, width = 7)
}

# plot_go(pheno = pheno, 
#         norm = norm, 
#         genes = as.vector(mapping$mgi_symbol), 
#         title = "BMP signaling pathway")

# plot_go(pheno = pheno, 
#         norm = norm, 
#         genes = as.vector(mapping$mgi_symbol), 
#         title = "TGF beta receptor signaling pathway")

# plot_go(pheno = pheno, 
#         norm = norm, 
#         genes = as.vector(mapping$mgi_symbol), 
#         title = "Notch signaling pathway")

# plot_go(pheno = pheno, 
#         norm = norm, 
#         genes = as.vector(mapping$mgi_symbol), 
#         title = "Wnt signaling pathway")

# plot_go(pheno = pheno, 
#         norm = norm, 
#         genes = as.vector(mapping$mgi_symbol), 
#         title = "ERK1 and ERK2 cascade")

plot_go(pheno = pheno, 
        norm = norm, 
        genes = as.vector(mapping$mgi_symbol), 
        title = "MAP kinase activity")



