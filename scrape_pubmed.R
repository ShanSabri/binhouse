## Setup
rm(list=ls(all=TRUE))
options(warn=-1)
pacman::p_load(tidyverse, RISmed, reshape2)
setwd("~/.Trash")


## Data
authors <- c("Kathrin Plath", "Shan Sabri", "Justin Langerman", "William Lowry")


## Scrape
results <- lapply(authors, function(x){
  message(x)
  res <- EUtilsSummary(x, type = "esearch", db = "pubmed")
  auths <- do.call(rbind.data.frame, Author(EUtilsGet(res)))
  auths %>% 
    separate(ForeName, c("ForeName", "initial")) %>%
    mutate(Author = paste(ForeName, LastName)) %>% 
    group_by(Author) %>% 
    summarise(n = n()) %>% 
    arrange(desc(n)) %>% 
    filter(Author != x) %>% 
    mutate(Query = x) -> to_ret
  return(to_ret)
}); names(results) <- authors


## Aggregrate 
tbl <- dcast(Author ~ Query, value.var = "n", data = do.call(rbind.data.frame, results))
tbl[is.na(tbl)] = 0 
row.names(tbl) <- tbl$Author; tbl$Author <- NULL


## Output
write.table(tbl, file.path(getwd(), "results.txt"), col.names = NA, sep = "\t", quote = FALSE)
