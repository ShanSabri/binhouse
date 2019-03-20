### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/other/MEF-ESC-program-boxplots/")




### LIBS
pacman::p_load(ggplot2, reshape2)
theme_set(theme_bw(base_size = 14) +
            theme(panel.grid = element_blank(),
                  legend.justification = "top"))




### BULK DATA
bulk <- read.table("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/kostas_rna/GSE90894-RPKM.txt", 
                   sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
row.names(bulk) <- toupper(row.names(bulk))
norm <- log2(bulk + 1)
norm <- norm[,c("MEFs", "ESCs")]
names(norm) <- gsub("s", "", names(norm))




### PROGRAMS
prog <- readRDS("~/Dropbox/PlathLab/Analyses/Justin_Langerman/Timecourse_analysis/data/metanetwork/2019.02.08.Drop10.All.Program.List.Cleanedx4.rds")
programs_set <- split(prog, f = prog$Set)
programs_run <- split(prog, f = prog$Run)
programs <- c(programs_run, programs_set)
lapply(programs, nrow)



### BOXPLOTS
facet_data <- lapply(sort(names(programs)), function(x){
  
  # x = "CDH1_P_8"
  
  program_genes <- as.vector(programs[[x]]$Genes)
  common_genes <- intersect(program_genes, row.names(norm))
  
  mtx <- norm[common_genes, , drop = FALSE]
  mtx$RATIO <- mtx$ESC - mtx$MEF
  
  to_plot <- reshape2::melt(mtx)
  
  ggplot(to_plot, aes(x = variable, y = value, fill = variable)) +
    geom_boxplot() +
    labs(x = "", y = "Log2(RPKM + 1)", 
         title = x) +
    theme(legend.position = "none") +
    scale_fill_manual(values=c("blue3", "red3", "green3")) +
    # ggsave(file.path(getwd(), "individual-prog", paste0("BOXPLOT-", x, ".png")), height = 6, width = 4)
    ggsave(file.path(getwd(), "individual-prog-with-ratio", paste0("BOXPLOT-", x, ".png")), height = 6, width = 5)
  
  to_plot$PROG <- x
  return(to_plot)
})

# FACET
facet_data <- do.call(rbind.data.frame, facet_data)
lvls <- gtools::mixedsort(unique(facet_data$PROG))
facet_data$PROG <- factor(facet_data$PROG, levels = lvls)
ggplot(facet_data, aes(x = variable, y = value, fill = variable)) +
  facet_wrap(~ PROG) +
  geom_boxplot() +
  labs(x = "", y = "Log2(RPKM + 1)", 
       title = "All programs") +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("blue3", "red3", "green3")) +
  # ggsave(file.path(getwd(), paste0("BOXPLOTS-FACET.png")), height = 10, width = 12)
  ggsave(file.path(getwd(), paste0("BOXPLOTS-FACET-WITH-RATIO.png")), height = 10, width = 15)

progs_to_show <- as.vector(unique(facet_data$PROG)[grepl("CDH1_|TC1_4.0", unique(facet_data$PROG))])
ggplot(subset(facet_data, PROG %in% progs_to_show), 
       aes(x = variable, y = value, fill = variable)) +
  facet_wrap(~ PROG, ncol = 5) +
  geom_boxplot() +
  labs(x = "", y = "Log2(RPKM + 1)") +
  theme(legend.position = "none") +
  scale_fill_manual(values=c("blue3", "red3", "green3")) +
  # ggsave(file.path(getwd(), paste0("BOXPLOTS-SUBSET-FACET.png")), height = 5, width = 10)
  ggsave(file.path(getwd(), paste0("BOXPLOTS-SUBSET-FACET-WITH-RATIO.png")), height = 5, width = 10)

# only ratio 
ratio_data <- subset(facet_data, variable == "RATIO")
lvls <- gtools::mixedsort(unique(ratio_data$PROG))
ratio_data$PROG <- factor(ratio_data$PROG, levels = lvls)

ggplot(ratio_data, aes(x = PROG, y = value)) +
  geom_boxplot(fill = "green3") +
  labs(x = "", y = "Log2(ESC_RPKM + 1) - Log2(MEF_RPKM + 1)", 
       title = "All programs", subtitle = "RATIO = ESC - MEF") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave(file.path(getwd(), paste0("BOXPLOTS-FACET-ONLY-RATIO.png")), height = 5, width = 15)

progs_to_show <- as.vector(unique(ratio_data$PROG)[grepl("CDH1_|TC1_4.0", unique(ratio_data$PROG))])
ggplot(subset(ratio_data, PROG %in% progs_to_show), aes(x = PROG, y = value)) +
  geom_boxplot(fill = "green3") +
  labs(x = "", y = "Log2(ESC_RPKM + 1) - Log2(MEF_RPKM + 1)",
       title = "Select programs", subtitle = "RATIO = ESC - MEF") +
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggsave(file.path(getwd(), paste0("BOXPLOTS-SUBSET-FACET-ONLY-RATIO.png")), height = 5, width = 8)




