### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/cpg-geneset-overlap-with-programs/")



### LIBS
pacman::p_load(ggplot2, reshape2)
theme_set(theme_bw(base_size = 14) + 
            theme(panel.grid = element_blank(), 
                  legend.justification = "top"))



### DATA
# PROGRAMS
data_dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data"
metanetwork <- readRDS(file.path(data_dir, "metanetwork", "2019.01.23.Drop10.All.Program.List.Cleanedx2.v2.rds"))
metanetwork <- split(metanetwork, f = metanetwork$Set)


# GENESETS
geneset_dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/Clustering/gene-lists/"
hcp <- toupper(as.vector(read.table(file.path(geneset_dir, "HCP_genename.txt"), header = FALSE)$V1))
icp <- toupper(as.vector(read.table(file.path(geneset_dir, "ICP_genename.txt"), header = FALSE)$V1))
lcp <- toupper(as.vector(read.table(file.path(geneset_dir, "LCP_genename.txt"), header = FALSE)$V1))



### COMPUTE OVERLAP
overlap <- do.call(rbind.data.frame, lapply(seq_along(metanetwork), function(x){
  id <- names(metanetwork[x]); message(id)
  prog <- as.vector(metanetwork[[x]]$Genes)
  prog_overlap <- data.frame(prog = id, 
                             hcp_overlap = (length(intersect(hcp, prog)) / length(prog)) * 100, 
                             icp_overlap = (length(intersect(icp, prog)) / length(prog)) * 100, 
                             lcp_overlap = (length(intersect(lcp, prog)) / length(prog)) * 100)
}))
overlap_m <- reshape2::melt(overlap, id.vars = "prog")



### PLOT
head(overlap_m)
ggplot(overlap_m, aes(x = prog, y = value, fill = variable)) + 
  facet_wrap(~variable, ncol = 1) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Percent Overlap") +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9)) +
  ggsave(file.path(getwd(), "PROG-PCT-OVERLAP.png"), height = 7, width = 9)
  
