## OPTS
rm(list=ls(all=TRUE))
setwd("~/Dropbox/PlathLab/Analyses/Lydia_Lee/2019-02-07/"); getwd()



## LIBS
pacman::p_load(Seurat, reshape2, 
               tidyverse, ggplot2, 
               reticulate)
theme_set(theme_bw(base_size=16))



## FXNS
filter_genes <- function(x, n_cells=5, is_expr=1) {
  return(x[rowSums(x >= is_expr) >= n_cells,])
}



## DATA
dges <- list.files(file.path(getwd(), "from_nathan_phillips", "delivered_1_31_19", "Expression Matrix"), 
                   pattern = "*.csv.gz", full.names = TRUE, recursive = TRUE)



# FILTER AND RUN SEURAT
objs <- lapply(dges, function(x){
  
  # data
  id <- sapply(strsplit(basename(x), "_", fixed=TRUE), function(x) (x[1]))
  out <- file.path(getwd(), "output", id)
  dir.create(out, showWarnings = FALSE)
  counts <- read.table(gzfile(x), header = TRUE, sep = ",", row.names = 1, check.names = FALSE)
  names(counts) <- paste(id, names(counts), sep = "_")
  counts <- filter_genes(counts)
  message(paste0(id, ": ", ncol(counts), " cells"))
  
  pheno <- data.frame(row.names = names(counts), 
                      ID = rep(id, length(names(counts))), 
                      Genes = apply(counts, 2, function(x) sum(x >= 1)),
                      UMIs = colSums(counts))
  
  # filter  
  before_n_cells <- nrow(pheno)
  ggplot(melt(pheno, id.vars = "ID"), aes(value)) +
    facet_wrap(~ variable, nrow = 1, scales = "free") +
    geom_density(aes(y = ..scaled.., fill = variable)) +
    labs(x = "", y = "Density", title = paste0("ID ", id, ", ", before_n_cells, " cells")) + 
    theme_bw(base_size = 14) +
    theme(legend.position = "none", 
          panel.grid = element_blank(), 
          strip.text = element_text(size = 18)) +
    ggsave(file.path(out, paste0(id, "-DIST-GENES-UMIS.png")), height = 4, width = 10)
  
  pheno <- subset(pheno, Genes >= 200)
  pheno <- subset(pheno, UMIs >= 500 & UMIs <= 30000)
  after_n_cells <- nrow(pheno)
  ggplot(melt(pheno, id.vars = "ID"), aes(value)) +
    facet_wrap(~ variable, nrow = 1, scales = "free") +
    geom_density(aes(y = ..scaled.., fill = variable)) +
    labs(x = "", y = "Density", title = paste0("ID ", id, ", ", after_n_cells, " cells (filtered)")) + 
    theme_bw(base_size = 14) +
    theme(legend.position = "none", 
          panel.grid = element_blank(), 
          strip.text = element_text(size = 18)) +
    ggsave(file.path(out, paste0(id, "-DIST-GENES-UMIS-FILT.png")), height = 4, width = 10)
  
  counts <- counts[, row.names(pheno)]
  counts <- counts[rowSums(counts) > 0,]
  identical(row.names(pheno), colnames(counts))  # sanity check
  
  # seurat 
  obj <- CreateSeuratObject(raw.data = counts, min.cells = 0, min.genes = 0, project = id)
  obj <- AddMetaData(object=obj, metadata=pheno)
  obj <- NormalizeData(object=obj, normalization.method="LogNormalize", scale.factor=10000)
  obj <- FindVariableGenes(object=obj, mean.function=ExpMean, dispersion.function=LogVMR, x.low.cutoff=0.0125, x.high.cutoff=3, y.cutoff=0.5); length(obj@var.genes)
  obj <- ScaleData(object=obj, vars.to.regress=c("nUMI"))
  # obj <- RunPCA(object=obj, pc.genes=obj@var.genes, do.print=FALSE)
  obj <- RunPCA(object=obj, pc.genes=row.names(obj@data), do.print=FALSE)
  obj <- FindClusters(object=obj, reduction.type="pca", dims.use=1:10, resolution=seq(0.2, 1, by=0.2), print.output=0, save.SNN=TRUE, temp.file.location=tempdir())
  obj <- RunTSNE(object=obj, dims.use=1:10, do.fast=FALSE)
  # obj <- RunTSNE(object=obj, dims.use=1:10, do.fast=TRUE, reduction.name="tsne3d", dim.embed=3)
  saveRDS(obj, file=file.path(out, paste0(id, "_seurat_obj.rds")), compress=TRUE)
  
  return(obj)
})
saveRDS(objs, file=file.path(getwd(), "output", "LIST-OF-SEURAT-OBJS.rds"), compress=TRUE)



## PLOT CLUSTERS
objs <- readRDS(file.path(getwd(), "output", "LIST-OF-SEURAT-OBJS.rds"))
lapply(seq_along(objs), function(x){
  out <- file.path(getwd(), "output", x)
  sapply(seq(0.2, 1, by=0.2), function(i){
    grp <- paste0("res.",i)
    message(paste(x, grp))
    tsne_coords <- GetDimReduction(objs[[x]], reduction.type = "tsne", slot = "cell.embeddings")
    tmp <- cbind(objs[[x]]@meta.data, tsne_coords)
    labs <- tmp %>%
      group_by_at(vars(one_of(grp))) %>%
      summarise(mean_x = mean(tSNE_1), mean_y = mean(tSNE_2))
    colnames(labs)[1] <- 'C'
    g <- ggplot(tmp, aes(tSNE_1, tSNE_2)) +
      geom_point(aes_string(colour=grp), size=0.5) +
      geom_label(data=labs, aes(mean_x, mean_y, label=C), hjust=0.5, vjust=0.5, color='black', size=8, label.padding=unit(0.10, "lines")) +
      labs(title=paste0('Resolution: ', i)) +
      theme(legend.position = 'none', panel.grid = element_blank()) +
      ggsave(file.path(out, paste0(unique(tmp$ID), '-TSNE-CLUSTER-RES-', i ,'.png')), height=6, width=7)
  })
})



## UMAP
use_python("/Users/shansabri/miniconda3/bin/python")
objs <- readRDS(file.path(getwd(), "output", "LIST-OF-SEURAT-OBJS.rds"))
lapply(seq_along(objs), function(y){
  # y = 1; x = 20; n = 15
  out <- file.path(getwd(), "output", y)
  obj <- objs[[y]]
  obj@var.genes <- row.names(obj@raw.data)  
  lapply(seq(10, 20, 5), function(x){
    lapply(seq(15, 25, 5), function(n){
      message(paste0("ID ", y, " -- PARAMS ", paste(x, n)))
      obj <- RunUMAP(obj, dims.use = 1:x, metric = "euclidean", seed.use = 2020, reduction.name = 'UMAP', n_neighbors = n)
      tmp <- cbind(obj@meta.data, GetDimReduction(obj, reduction.type = 'UMAP', slot = "cell.embeddings"))
      # saveRDS(obj, file=file.path(out, paste0(y, "_seurat_obj.rds")), compress=TRUE)
      
      # sapply(seq(0.2, 1, by=0.2), function(i){
        # grp <- paste0("res.",i)
        # message(grp)
        # labs <- tmp %>%
          # group_by_at(vars(one_of(grp))) %>%
          # summarise(mean_x = mean(UMAP1), mean_y = mean(UMAP2))
        # colnames(labs)[1] <- 'C'
        g <- ggplot(tmp, aes(UMAP1, UMAP2)) +
          # geom_point(aes_string(colour=grp), size=0.5) +
          geom_point(size=0.5, colour = "blueviolet") +
          # geom_label(data=labs, aes(mean_x, mean_y, label=C), hjust=0.5, vjust=0.5, color='black', size=8, label.padding=unit(0.10, "lines")) +
          # labs(title=paste0('Resolution: ', i)) +
          theme(legend.position = 'none', panel.grid = element_blank()) +
          # ggsave(file.path(out, paste0(y, "-UMAP-TC1-TIMEPOINT-",x,"-NEIGHBORS-",n,"-RES-",i,".png")), height=6, width=7)
          ggsave(file.path(out, paste0(y, "-UMAP-TC1-TIMEPOINT-",x,"-NEIGHBORS-",n,".png")), height=6, width=7)
      # })
    })
  })
})



# ## DIFF EXP
# obj <- readRDS(paste(getwd(), paste0(ID, '_seurat_obj.rds'), sep='/'))
# idx <- grep("res.", colnames(obj@meta.data))
# res <- sort(colnames(obj@meta.data)[idx])
# 
# mclapply(res, function(r){
#   message(r)
#   obj  <- SetAllIdent(obj, id=r)
#   message('...running FindAllMarkers()')
#   tmp  <- FindAllMarkers(obj, test.use="bimod", only.pos=TRUE, min.pct=0.25, thresh.use=0.25)
#   outf <- paste(getwd(), paste0(ID, '_diffexp_',r,'_bimod.rds'), sep='/')
#   message(outf)
#   saveRDS(tmp, file=outf)
# }, mc.cores = 7)