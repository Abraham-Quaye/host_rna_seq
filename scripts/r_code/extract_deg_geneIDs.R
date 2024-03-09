#!/usr/bin/env Rscript

library(tidyverse)

source("scripts/r_code/deg_analysis.R")

# write.table(bg_gene_ids, file = "results/r/tables/all_bg_geneIDs.txt",
#             quote = F, row.names = F, col.names = F)

# function to filter out DEG regulation types at different timepoints
xtract_geneID_tp <- function(data = diff_exp_genes, reg, tp){
  ids <- data %>%
    filter(timepoint == tp, regulation == reg) %>%
    distinct(gene_id)
  
  return(ids)
}

# Function to extract all up/down-regulated genes
xtract_all_deg_regTypes <- function(data = diff_exp_genes, reg_type){
  ids <- data %>%
  filter(regulation == reg_type) %>%
  distinct(gene_id)
  return(ids)
}

all_deg_tables <- list(
  # Extract total background genes
  all_bg_geneIDs = crude_diff_exp_genes %>%
    distinct(gene_id),
  # Extract 4hr down-regulated genes
  down_degIDs_4hrs = xtract_geneID_tp(reg = "down", tp = "t4"),
  # Extract 12hr up/down-regulated genes
  down_degIDs_12hrs = xtract_geneID_tp(reg = "down", tp = "t12"),
  up_degIDs_12hrs = xtract_geneID_tp(reg = "up", tp = "t12"),
  # Extract 24hr up/down-regulated genes
  down_degIDs_24hrs = xtract_geneID_tp(reg = "down", tp = "t24"),
  up_degIDs_24hrs = xtract_geneID_tp(reg = "up", tp = "t24"),
  # Extract all up/down-regulated genes
  all_down_degs = xtract_all_deg_regTypes(reg_type = "down"),
  all_up_degs = xtract_all_deg_regTypes(reg_type = "up")
)

# Save geneID tables
for(t in names(all_deg_tables)){
  print(paste("Saving:", t))
  write.table(all_deg_tables[[t]],
              file = paste0("results/r/tables/", t, ".txt"),
              quote = F, row.names = F, col.names = F)
}