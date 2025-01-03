#!/usr/bin/env Rscript

library(tidyverse)


# # Download GO data from DAVID online analysis resource
# get_keggData <- function(link, filename){
#   download.file(url = link,
#                 destfile = paste0("results/r/tables/davidGO_", filename, ".tsv"))
# }
# 
#=================== 12hpi data
down12bp <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491024564.txt"
down12cc <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491078915.txt"
down12mf <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491113661.txt"

up12bp <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491266618.txt"
up12cc <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491292554.txt"
up12mf <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491329826.txt"

#=================== 24hpi data
down24bp <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491482223.txt"
down24cc <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491561991.txt"
down24mf <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491589433.txt"

up24bp <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491686088.txt"
up24cc <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491885943.txt"
up24mf <- "https://david.ncifcrf.gov/data/download/chart_26C68EFBFD081725491936175.txt"

go_data <- list(down12bp = down12bp, down12cc = down12cc, down12mf = down12mf,
                up12bp = up12bp, up12cc = up12cc, up12mf = up12mf,
                down24bp = down24bp, down24cc = down24cc, down24mf = down24mf,
                up24bp = up24bp, up24cc = up24cc, up24mf = up24mf)
# 
# map2(go_data, as.list(names(go_data)), ~get_keggData(link = .x, filename = .y))


###############################################################################
#             Process GO results from DAVID to include as table in paper

process_DAVID_GOresults <- function(go_result){
  pdata <- read_tsv(go_result, col_names = T,
                    col_select = c(source = Category, term = Term, num_degs = Count,
                                   fold_enrich = `Fold Enrichment`, recall = `%`,
                                   pval = PValue, qval = Benjamini)) %>%
    mutate(term = str_replace(term, "^mgp\\d+:([a-zA-Z\\s-\\(\\)/]+)$", "\\1"),
           recall = recall / 100) %>%
    filter(qval <= 0.05) %>%
    mutate(across(where(is.numeric), ~format(.x, digits = 3)))
  
  return(pdata)
}

david_go_files <- list.files(path = "results/r/tables", full.names = T,
                             pattern = "davidGO_(up|down)\\d{2}(bp|cc|mf)\\.tsv") %>%
  set_names(sort(names(go_data)))

david_go_results <- map_dfr(david_go_files, process_DAVID_GOresults, .id = "timepoint") %>%
  mutate(term = str_replace(term, "GO:\\d+~(.+)", "\\1"),
         across(c(num_degs:qval), as.numeric),
         source = case_match(source,
                             "GOTERM_BP_ALL" ~ "GO:BP",
                             "GOTERM_CC_ALL" ~ "GO:CC",
                             "GOTERM_MF_ALL" ~ "GO:MF")) %>%
  filter(qval <= 0.05)


prep_GO_table <- function(go_res){
  # Ensure GOTerms are ordered as: BP, CC, and MF
  go_res <- go_res %>% arrange(source)
  
  # Calculate row number to insert GOTerm label
  GObp_row_num <-  go_res %>% filter(source == "GO:BP") %>% nrow(.)
  GOcc_row_num <- go_res %>% filter(source == "GO:CC") %>% nrow(.) + GObp_row_num
  
  # GOTerm labels to insert into table
  bp <- " Biological Process"
  cc <- " Cellular Component"
  mf <- " Molecular Function"
  
  # Insert GOTerm labels into table and return results
  res_table <- go_res %>%
    add_row(timepoint = bp, pval = NA, num_degs = 0, recall = 0,
            source = "GO:BP", term = bp, .before = 1) %>%
    add_row(timepoint = cc, pval = NA, num_degs = 0, recall = 0,
            source = "GO:CC", term = cc, .after = GObp_row_num + 1) %>% 
    add_row(timepoint = mf, pval = NA, num_degs = 0, recall = 0,
            source = "GO:MF", term = mf, .after = GOcc_row_num + 2) %>%
    select(source, term, fold_enrich, num_degs, qval) %>%
    mutate(qval = base::format(qval, scientific = T, digits = 3))
  
  return(res_table)
}


tab_res <- tibble(timepoint = c("t12", "t24"),
                  res = list(filter(david_go_results, str_detect(timepoint, "12")),
                             filter(david_go_results, str_detect(timepoint, "24"))),
                  down = map(res, ~ filter(.x, str_detect(timepoint, "^down"))),
                  up = map(res, ~ filter(.x, str_detect(timepoint, "^up"))),
                  go_tab_down = map(down, prep_GO_table),
                  go_tab_up = map(up, prep_GO_table))

## KEGG RESULTS FROM Gprofiler2 ================

read_enrichment_data <- function(go_res_file, type){
  res <- read_tsv(go_res_file, col_names = T,
           col_select = -c(significant, term_size, term_id,
                           query_size, precision, recall,
                           effective_domain_size:parents)) %>%
    filter(str_detect(source, type)) %>%
    split(.$query)
  
  return(res)
}

t12_kegg_res <- read_enrichment_data("results/r/tables/t12_GO_results.tsv", "KEGG")
t24_kegg_res <- read_enrichment_data("results/r/tables/t24_GO_results.tsv", "KEGG")


kegg_res_all <- bind_rows(c(t12_kegg_res, t24_kegg_res))