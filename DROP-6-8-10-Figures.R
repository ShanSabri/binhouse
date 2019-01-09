### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/Drop-6-8-10/")



### LIBS
pacman::p_load(ggplot2, reshape2, RColorBrewer)
theme_set(theme_bw(base_size = 14) + 
            theme(panel.grid = element_blank(), 
                  legend.justification = "top"))
pal1 <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")




### FXNS
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  substr(x, 2, nchar(x)) <- tolower(substr(x, 2, nchar(x)))
  return(x)
}

range_norm <- function(x){
  (x-min(x)) / (max(x)-min(x))
}




### DATA 
data_dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data"
raw <- readRDS(paste(data_dir, "DROP-6-8-10-TOMATO-RAW-FILTERED-DGE.rds", sep = "/"))
norm <- log((sweep(raw, 2, apply(raw, 2, function(x) sum(x)), "/")*10000)+1)
pheno <- readRDS(paste(data_dir, "DROP-6-8-10-TOMATO-PHENODATA.rds", sep = "/"))
identical(names(norm), row.names(pheno))



### PROGRAMS
programs <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/Drop10.Aggregated.Programs.Gene.Master.List.rds")
programs <- split(programs, f = programs$Set)



### PLOT PROGRAMS
lapply(seq_along(programs), function(x){
  id <- names(programs)[x]
  p <- programs[[x]]
  message(id)
  
  md <- t(subset(norm[, row.names(pheno)], row.names(norm) %in% as.vector(p$Genes)))
  md <- cbind.data.frame(pheno, exp = range_norm(rowMeans(md)))  
  md$Run <- factor(md$Run, levels = c("D6", "D8", "D10"))
  
  ggplot(md[with(md, order(exp)), ], aes(tSNE_1, tSNE_2)) +
    facet_wrap(~Run, nrow = 1) +
    geom_point(data = md[,c("tSNE_1", "tSNE_2")], aes(tSNE_1, tSNE_2), colour = "grey") +
    geom_point(aes(colour = exp)) +
    scale_colour_gradientn(colours = c("grey70","grey50","lightgoldenrod","orange","red"), name = "") +
    labs(title = id) +
    theme(
      axis.title = element_text(size = 28),
      axis.text = element_text(size = 16),
      plot.title = element_text(hjust = 0.5, size = 32),
      plot.subtitle = element_text(hjust = 0.5, size = 20),
      legend.text = element_text(size = 24),
      legend.title = element_blank(),
      legend.position = "top",
      legend.key.width = unit(2, "cm"),
      strip.text = element_text(size = 22)
    ) +
    guides(colour = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
    # ggsave(paste(getwd(), "programs-avg-exp", paste0("TSNE-PROGRAM-",id,".png"), sep = "/"), height = 8, width = 7)
    ggsave(paste(getwd(), "programs-avg-exp", "facet-run", paste0("TSNE-PROGRAM-",id,"-FACET.png"), sep = "/"), height = 6, width = 11)
})


### PLOT TSNE MARKER GENES
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
    ggsave(paste(getwd(), "markers", paste0("TSNE-", g, ".png"), sep = "/"), height = 8, width = 7)
  return(p)
}
genes <- c("AVIL", "BCL11A", "CACNA2D2", "CRYM", "CHST3", "EHF", "NT5E", "ITGB4", "ITGA4")
lapply(genes, plot_single_gene)



### TSNE 
p <- pheno
names(p)[names(p)=='Timepoint'] = 'Day'
p$Day <- gsub('Day ', '', p$Day)
table(p$Day)
p$Day <- factor(p$Day, levels=c(0,3,6,7,9,12,13,15))
ggplot(p, aes(tSNE_1, tSNE_2)) +
  geom_point(aes(colour=Day)) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size=4))) + 
  ggsave(paste(getwd(), "TSNE.png", sep="/"), height=6, width=6.5)




### TSNE - FACET
ggplot(pheno, aes(tSNE_1, tSNE_2)) +
  facet_wrap(~Timepoint, nrow = 2) +
  geom_point(data=pheno[,c('tSNE_1', 'tSNE_2')], colour="grey", size=0.5) +  
  geom_point(aes(colour=Timepoint), size=0.5) +
  theme(legend.position="none") +
  theme(strip.text.x = element_text(size = 18)) +
  theme(axis.title=element_text(size=14), axis.text = element_text(size=10)) +
  ggsave(paste(getwd(), "TSNE-FACET.png", sep="/"), height=6, width=11)




### TSNE - COL1A1, NANOG, ESRRB
genes_of_interest <- c("NANOG", "ESRRB", "COL1A1")
marker_data <- cbind(pheno[,c('Timepoint','tSNE_1','tSNE_2')], t(subset(norm[,row.names(pheno)], row.names(norm) %in% toupper(genes_of_interest))))
marker_data_m <- reshape2::melt(marker_data, id.vars = c("Timepoint", "tSNE_1", "tSNE_2"))
marker_data_m <- rbind(
  marker_data_m,
  data.frame(Timepoint=rep('All', nrow(marker_data)), tSNE_1=marker_data$tSNE_1, tSNE_2=marker_data$tSNE_2, variable=rep('NANOG', nrow(marker_data)), value=marker_data$NANOG),
  data.frame(Timepoint=rep('All', nrow(marker_data)), tSNE_1=marker_data$tSNE_1, tSNE_2=marker_data$tSNE_2, variable=rep('ESRRB', nrow(marker_data)), value=marker_data$ESRRB),
  data.frame(Timepoint=rep('All', nrow(marker_data)), tSNE_1=marker_data$tSNE_1, tSNE_2=marker_data$tSNE_2, variable=rep('COL1A1', nrow(marker_data)), value=marker_data$COL1A1)
)
marker_data_m$variable <- factor(marker_data_m$variable, levels=c(genes_of_interest, 'All'))
p <- marker_data_m[with(marker_data_m,order(variable, value)),]
# p$value[p$value > 3.5] <- 3.5
ggplot(p, aes(tSNE_1, tSNE_2, colour=value)) +
  facet_grid(variable~Timepoint, switch = 'y') +
  geom_point(data=marker_data[,c("tSNE_1", "tSNE_2")], aes(tSNE_1, tSNE_2), colour="grey", size=1) +  
  geom_point() +
  scale_color_gradientn(colours = pal1(100), 'Normalized\nExpression') + 
  scale_y_continuous(position = "right") +
  theme(strip.text.x = element_text(size = 24), 
        strip.text.y = element_text(size=24),
        axis.title=element_text(size=20),
        axis.text = element_text(size=14), 
        legend.title=element_text(size=22),
        legend.position="top", 
        legend.justification="left", 
        legend.direction="horizontal", 
        text = element_text(size=12),
        legend.key.size=unit(1,"cm"), 
        legend.text = element_text(size=16)) +
  ggsave(paste(getwd(), "TSNE-COL1A1-NANOG-ESRRB.png", sep="/"), height=7, width=22)


