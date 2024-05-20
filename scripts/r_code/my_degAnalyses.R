#!/usr/bin/env Rscript

library(tidyverse)
library(readxl)
library(magrittr)

# load in files for differentially expressed genes
diff_files <- list.files("results/r/tables",
                         pattern = "signif_\\d{1,2}hrsDEGs\\.csv",
                         full.names = T) %>%
  setNames(c("t12", "t24", "t4", "t72"))

# function to read in all the files in a compatible manner
pull_exp_data <- function(csv){
  df <- read.csv(file = csv)
  
  # insert column of NA values for the third infected replicate
  if(csv == diff_files["t12"]){
    df <- df %>%
      add_column("12" = NA, .after = 2)
  }
  
  # set the column names
  df <- df %>%
    set_colnames(c("gene_id", "inf_fpkm_r1", "inf_fpkm_r2", "inf_fpkm_r3", "mock_fpkm_r1",
                   "mock_fpkm_r2", "log2fc", "pval", "qval"))
  
  return(df)
}

# load all differentially expressed genes from all time points
crude_diff_exp_genes <- map_dfr(diff_files, pull_exp_data, .id = "timepoint")

diff_exp_genes <- crude_diff_exp_genes %>%
  filter(qval <= 0.05) %>%
  mutate(regulation = ifelse(log2fc > 0, "up", "down"))
