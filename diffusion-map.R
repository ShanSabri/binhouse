### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/diffusion-map")




### LIBS
pacman::p_load(ggplot2, reshape2, biomaRt)
theme_set(theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.justification = "top"
  ))




### DATA
data_dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data"
norm <- readRDS(paste(data_dir, "TC1-NORM-FILTERED-DGE.rds", sep = "/"))
pheno <- readRDS(paste(data_dir, "TC1-PHENODATA.rds", sep = "/"))
identical(names(norm), row.names(pheno))




### SCALE + PCA
new_scale <- function(x, do_scale = TRUE, do_center = TRUE) {
  bin_size <- 1000
  max_bin <- floor(nrow(x) / bin_size) + 1
  for (i in 1:max_bin) {
    my_inds <- ((bin_size * (i - 1)):(bin_size * i - 1)) + 1
    my_inds <- my_inds[my_inds <= nrow(x)]
    new_data <- t(scale(t(as.matrix(x[my_inds, ])), center = do_center, scale = do_scale))
    new_data[new_data > 10] <- 10
    x[my_inds, ] <- new_data
  }
  x
}
scaled <- new_scale(norm)
pca <- irlba::irlba(t(scaled), nv = 30)
projected_loadings <- data.frame(t(scaled) %*% as.matrix(pca$v))
projected_loadings$Timepoint <- pheno$Timepoint
# identical(row.names(projected_loadings), row.names(pheno))




### DIFFUSION MAP
n_k <- 50
dm <- destiny::DiffusionMap(projected_loadings,
  k = n_k,
  n_eigs = 3,
  verbose = TRUE
)
saveRDS(dm, compress = TRUE, file.path(getwd(), paste0("DM-OBJ-K", n_k, ".rds")))


get_louvain_clusters <- function(transitions) {
  graph <- igraph::graph_from_adjacency_matrix(transitions, 'undirected', weighted = TRUE)
  as.integer(unclass(igraph::membership(igraph::cluster_louvain(graph))))
}
DM_cluster_assignments <- factor(get_louvain_clusters(dm@transitions))



### PLOT
gg_color_hue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

rename <- function(x) {
  if (x < 10) {
    return(name <- paste("000", i, "plot.png", sep = ""))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste("00", i, "plot.png", sep = ""))
  }
  if (x >= 100) {
    return(name <- paste("0", i, "plot.png", sep = ""))
  }
}

cols <- data.frame(
  DAY = levels(pheno$Timepoint),
  COLOR = gg_color_hue(length(levels(pheno$Timepoint)))
)
to_plot <- cbind.data.frame(pheno, dm@eigenvectors, DM_cluster_assignments)



# 2D
A <- ggplot(to_plot, aes(DC1, DC2)) + geom_point(aes(colour = Timepoint))
B <- ggplot(to_plot, aes(DC1, DC3)) + geom_point(aes(colour = Timepoint))
C <- ggplot(to_plot, aes(DC1, DC2)) + geom_point(aes(colour = DM_cluster_assignments)) + labs(colour = "Cluster") + scale_colour_viridis_d()
D <- ggplot(to_plot, aes(DC1, DC3)) + geom_point(aes(colour = DM_cluster_assignments)) + labs(colour = "Cluster") + scale_colour_viridis_d()
cowplot::plot_grid(A, B, C, D, labels = "AUTO") +
  ggsave(file.path(getwd(), paste0("DM-K", n_k, ".png")), height = 13, width = 15)
dev.off()



# 3D
to_plot$Timepoint_Colour <- cols$COLOR[match(to_plot$Timepoint, cols$DAY)]
png(file.path(getwd(), paste0("DM-K", n_k, "-3D.png")), height = 700, width = 900)
with(to_plot, scatterplot3d::scatterplot3d(
  x = DC1, y = DC2, z = DC3,
  color = Timepoint_Colour,
  pch = 19
))
legend("top",
  inset = -0.05, bty = "n", cex = 1, title = "", as.vector(unique(to_plot$Timepoint)),
  fill = as.vector(unique(to_plot$Timepoint_Colour)), xpd = TRUE, horiz = TRUE
)
dev.off()




# 3D ANIMATED
for(i in seq(1, 360, by = 2)){
  message(i)
  png(file.path(getwd(), rename(i)), height = 700, width = 900)
  with(to_plot, scatterplot3d::scatterplot3d(
    x = DC1,
    y = DC2,
    z = DC3,
    color = Timepoint_Colour,
    pch = 19,
    angle = i
  ))
  legend("top",
    inset = -0.05, bty = "n", cex = 1, title = "", as.vector(unique(to_plot$Timepoint)),
    fill = as.vector(unique(to_plot$Timepoint_Colour)), xpd = TRUE, horiz = TRUE
  )
  dev.off()
}
my_command <- paste("convert *plot.png -delay 1 -loop 0 ", file.path(getwd(), paste0("DM-K", n_k, "-3D-ANIMATED.gif")))
system(my_command)
cleanup <- "rm *plot*.png"
system(cleanup)

# 3D INTERACTIVE
# library(rgl)
# plot3d(destiny::eigenvectors(dm)[, 1:3], type = 's', radius = .01)
# view3d(theta = 10, phi = 30, zoom = .8)
# rgl.close()