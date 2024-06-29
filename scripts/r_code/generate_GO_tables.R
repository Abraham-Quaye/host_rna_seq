#!/usr/bin/env Rscript

library(tidyverse)

read_enrichment_data <- function(go_res_file, type){
  res <- read_tsv(go_res_file, col_names = T,
           col_select = -c(significant, term_size, term_id,
                           query_size, precision, recall,
                           effective_domain_size:parents)) %>%
    filter(str_detect(source, type)) %>%
    split(.$query)
  
  return(res)
}

t12_go_res <- read_enrichment_data("results/r/tables/t12_GO_results.tsv", "GO:")
t24_go_res <- read_enrichment_data("results/r/tables/t24_GO_results.tsv", "GO:")

t12_kegg_res <- read_enrichment_data("results/r/tables/t12_GO_results.tsv", "KEGG")
t24_kegg_res <- read_enrichment_data("results/r/tables/t24_GO_results.tsv", "KEGG")

prep_GO_table <- function(go_res){
  # Ensure GOTerms are ordered as: BP, CC, and MF
  go_res <- go_res %>% arrange(source)
  
  # Calculate row number to insert GOTerm label
  GObp_row_num <-  go_res %>% filter(source == "GO:BP") %>% nrow(.)
  GOcc_row_num <- go_res %>% filter(source == "GO:CC") %>% nrow(.) + GObp_row_num
  
  # GOTerm labels to insert into table
  bp <- "Biological Process"
  cc <- "Cellular Component"
  mf <- "Molecular Function"
  
  # Insert GOTerm labels into table and return results
  res_table <- go_res %>%
    add_row(query = bp, p_value = 0, intersection_size = 0, 
            source = bp, term_name = bp, .before = 1) %>%
    add_row(query = cc, p_value = 0, intersection_size = 0, 
            source = cc, term_name = cc, .after = GObp_row_num + 1) %>% 
    add_row(query = mf, p_value = 0, intersection_size = 0, 
            source = mf, term_name = mf, .after = GOcc_row_num + 2) %>%
    select(source, term_name, p_value, intersection_size) %>%
    mutate(p_value = base::format(p_value, scientific = T, digits = 3))
  
  return(res_table)
}

goTerm_tables <- list(t12_go_res_down = prep_GO_table(t12_go_res[["down_degIDs_12hrs"]]),
                      t12_go_res_up = prep_GO_table(t12_go_res[["up_degIDs_12hrs"]]),
                      t24_go_res_down = prep_GO_table(t24_go_res[["down_degIDs_24hrs"]]),
                      t24_go_res_up = prep_GO_table(t24_go_res[["up_degIDs_24hrs"]]),
                      kegg_res_all = bind_rows(c(t12_kegg_res, t24_kegg_res))
                      )
