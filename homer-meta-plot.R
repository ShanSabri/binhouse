### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/homer/TC1/homer_output/")




### LIBS
pacman::p_load(marge, tidyverse, 
               dplyr, gplots)




### AGGREGRATE RESULTS 
dir <- "~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/homer/TC1/homer_output/500_500"
res <- list.dirs(dir, full.names = TRUE, recursive = FALSE)
all_motif_results <- do.call(rbind.data.frame, lapply(seq_along(res), function(x) {
  progs <- list.dirs(res[[x]], full.names = TRUE, recursive = FALSE)
  if(length(progs) != 0){
    k <- lapply(progs, function(xx){
      known <- read_known_results(path = xx, homer_dir = TRUE)
      known$ID <- paste(basename(res)[x], basename(xx), sep = "_")
      # known <- subset(known, database == "Homer")
      known$motif_pwm <- NULL
      known
    })
  } 
  return(do.call(rbind.data.frame, k))
}))
all_motif_results <- split(all_motif_results, f = all_motif_results$ID)
saveRDS(all_motif_results, compress = TRUE, paste(dir, "RESULTS-AGGREGRATED.rds", sep = "/"))




### WHAT ARE THE TOP 10 MOTIFS FOR EACH GENE SET?
all_motif_results %>% 
  bind_rows() %>% 
  mutate(prog = sapply(strsplit(.$ID, "_", fixed=TRUE), function(x) (x[2])), 
         ID = sapply(strsplit(.$ID, "_", fixed=TRUE), function(x) (x[1]))) %>%  
  group_by(ID, prog) %>%
  select(ID, prog, motif_name, motif_family, experiment, accession, log_p_value) %>%
  top_n(log_p_value, n = 10) %>%
  ungroup() %>%
  write_tsv(paste(dir, "TOP10-MOTIFS-BY-PROGRAM.txt", sep = "/"), na = "NA",  quote_escape = "double")




### PLOT HEATMAP 
all_motif_results <- all_motif_results %>% bind_rows()
motif_df <- reshape2::dcast(all_motif_results, motif_name + motif_family + accession + consensus ~ ID,
                            value.var = "log_p_value", fun.aggregate = max)
row.names(motif_df) <- paste(motif_df$motif_name, motif_df$motif_family, motif_df$accession, motif_df$consensus, sep = "/")
motif_df$motif_name <- NULL
motif_df$motif_family <- NULL
motif_df$accession <- NULL
motif_df$consensus <- NULL

motif_df$max_p_val <- apply(motif_df, 1, max)

maxs <- seq(0, 8, by = 1)
# pdf(paste(getwd(), paste0("ENRICHED_MOTIFS_MAX_THRESHOLD.pdf"), sep = "/"), height = 20, width = 13)
pdf(paste(getwd(), paste0("ENRICHED_MOTIFS_MAX_THRESHOLD_LOG_LOG.pdf"), sep = "/"), height = 20, width = 13)
lapply(maxs, function(x) {
  # x = 0
  tmp <- subset(motif_df, max_p_val >= x)
  tmp$max_p_val <- NULL
  print(dim(tmp))
  heatmap.2(as.matrix(log10(tmp + 1)),
  # heatmap.2(as.matrix(tmp),
            Rowv = TRUE,
            Colv = FALSE,
            trace = "none",
            na.color = "grey",
            margins = c(20, 25),
            col = cm.colors(255),
            keysize=0.5,
            lhei = c(0.5, 6), 
            main = paste("Max log_p_value >=", x))
})
dev.off()