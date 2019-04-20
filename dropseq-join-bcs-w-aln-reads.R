### Setup
rm(list = ls(all = TRUE))
options(warn = -1)
setwd("~/Desktop/")



### Data
bcs_file <- file.path(getwd(), "AP35_bcs_nameSorted.sam")
bed_file <- file.path(getwd(), "AP35_hg19_Final_nameSorted.bed")
bcs <- tibble::as.tibble(data.table::fread(bcs_file))[, c("V1", "V10")]
bed <- tibble::as.tibble(data.table::fread(bed_file))



### Merge in barcodes 
output <- dplyr::left_join(x = bed, y = bcs, by = c("V4" = "V1"))


### Cleanup
rm(bcs, bed); gc()
output$BC <- paste0("BC:", substr(output$V10.y, 1, 12))
output$UMI <- paste0("UMI:", substr(output$V10.y, 13, 20))
output$V10.x <- paste0("GENE:", output$V10.x)
output$V11 <- paste0("ACC:", output$V11)
output$V1 <- output$V7 <- paste0("HUMAN_", output$V1)
output$FEAT <- "FEAT:CODING"
output <- output[, c("V4", "V1", "V2", "V3", "V6", "V7", "V8", "V9", "V10.x", 
                 "V11", "V12", "V10.y", "BC", "UMI", "FEAT")]


### Output
data.table::fwrite(output, file = file.path(getwd(), "AP35_hg19_Joined.bed"), 
                   col.names = FALSE, quote = FALSE, sep = "\t", showProgress = TRUE, verbose = TRUE)
