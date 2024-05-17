#!/usr/bin/env Rscript

library(tidyverse)
library(readxl)
library(magrittr)

# load in files for differentially expressed genes
diff_files <- list.files("raw_analysis/diff_exp",
                         pattern = "diff_gene_exp_\\d{1,2}hrs.xlsx",
                         full.names = T) %>%
  setNames(c("t12", "t24", "t4", "t72"))

# function to read in all the files in a compatible manner
pull_exp_data <- function(workbook){
  df <- read_excel(workbook, sheet = "Sheet1",
             col_names = T)
  
  # insert column of NA values for the third infected replicate
  if(workbook == diff_files["t12"]){
    df <- df %>%
      add_column("12" = NA, .after = 11)
  }
  
  # set the column names
  df <- df %>%
    set_colnames(c("gene_id", "gene_name", "transcript_id", "go_term",
                   "kegg", "ko_entry", "ec", "description", "trans_type",
                   "inf_fpkm_r1", "inf_fpkm_r2", "inf_fpkm_r3", "mock_fpkm_r1",
                   "mock_fpkm_r2", "fc", "log2fc", "pval", "qval", "regulation",
                   "signif"))
  
  return(df)
}

# load all differentially expressed genes from all time points
crude_diff_exp_genes <- map_dfr(diff_files, pull_exp_data, .id = "timepoint") %>%
  mutate(gene_id = ifelse(str_detect(gene_id, "gene-"),
                          str_replace(gene_id, "gene-LOC(\\d+)", "\\1"),
                          gene_id))

diff_exp_genes <- crude_diff_exp_genes %>%
  filter(qval <= 0.05) %>%
  dplyr::select(-signif)
