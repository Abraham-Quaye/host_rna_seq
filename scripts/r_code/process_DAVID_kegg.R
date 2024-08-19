#!/usr/bin/env Rscript

library(tidyverse)
library(magrittr)

# Download KEGG data from DAVID online analysis resource
# down12 <- "https://david.ncifcrf.gov/data/download/chart_3AC8819126C41724101963794.txt"
# up12 <- "https://david.ncifcrf.gov/data/download/chart_798F8A48EC0A1719613700031.txt"
# 
# down24 <- "https://david.ncifcrf.gov/data/download/chart_798F8A48EC0A1719614251153.txt"
# up24 <- "https://david.ncifcrf.gov/data/download/chart_798F8A48EC0A1719614284567.txt"
# 
# get_keggData <- function(link, filename){
#   download.file(url = link,
#                 destfile = paste0("results/r/tables/davidKEGG_", filename, "hrs.tsv"))
# }
# 
# map2(list(down12, up12, down24, up24),
#     list("down12", "up12", "down24", "up24"), ~get_keggData(link = .x, filename = .y))
# 

###############################################################################
#             Process KEGG results from DAVID to plot as table in paper

process_DAVID_kegg <- function(kegg_data){
  pdata <- read_tsv(kegg_data, col_names = T,
           col_select = c(term = Term, num_degs = Count, rich_factor = `%`,
                          fold_enrich = `Fold Enrichment`, pval = PValue,
                          qval = Benjamini)) %>%
    mutate(term = str_replace(term, "^mgp\\d+:([a-zA-Z\\s-\\(\\)/]+)$", "\\1")) %>%
    filter(qval <= 0.05) %>%
    mutate(across(where(is.numeric), ~format(.x, digits = 3)))
  
  return(pdata)
}
kegg_files <- list(down12hrs = "results/r/tables/davidKEGG_down12hrs.tsv",
                   up12hrs = "results/r/tables/davidKEGG_up12hrs.tsv",
                   down24hrs = "results/r/tables/davidKEGG_down24hrs.tsv",
                   up24hrs = "results/r/tables/davidKEGG_up24hrs.tsv")

kegg_result <- map_dfr(kegg_files, process_DAVID_kegg, .id = "timepoint")
