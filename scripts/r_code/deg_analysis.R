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

# # function to ensure consistency with groups/condition/treatment
check_grp_consistency <- function(grp_fpkms){
  # condition to reduce stringency (max(grp_fpkms, na.rm = T) > 1) & -- add to if-statement
  if(max(grp_fpkms, na.rm = T) >= (2*min(grp_fpkms, na.rm = T))){
    return(NA)
  }else{
    return(mean(grp_fpkms, na.rm = T))
  }
}

# # load all differentially expressed genes from all time points
crude_diff_exp_genes <- map_dfr(diff_files, pull_exp_data, .id = "timepoint") %>%
  mutate(gene_id = ifelse(str_detect(gene_id, "gene-"),
                          str_replace(gene_id, "gene-LOC(\\d+)", "\\1"),
                          gene_id))

diff_exp_genes <- crude_diff_exp_genes %>%
  mutate(inf_mean_fpkm = apply(.[, c("inf_fpkm_r1", "inf_fpkm_r2", "inf_fpkm_r3")], 1, check_grp_consistency),
         mock_mean_fpkm = apply(.[, c("mock_fpkm_r1", "mock_fpkm_r2")], 1, check_grp_consistency)) %>%
  drop_na(inf_mean_fpkm, mock_mean_fpkm) %>% 
# filter for DEGs
  filter(qval <= 0.05) %>%
  mutate(gene_known = ifelse(str_detect(gene_name, "LOC"), F, T))
