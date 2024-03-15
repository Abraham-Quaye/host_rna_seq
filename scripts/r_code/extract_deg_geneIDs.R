#!/usr/bin/env Rscript

source("scripts/r_code/deg_analysis.R")

# function to filter out DEG regulation types at different timepoints
xtract_geneID_tp <- function(data = diff_exp_genes, reg, tp){
  ids <- data %>%
    filter(timepoint == tp, regulation == reg) %>%
    select(gene_id, gene_name, log2fc) %>%
    as.data.frame() %>% 
    distinct(gene_id, gene_name, .keep_all = T) %>%
    set_rownames(.$gene_id)
  
  return(ids)
}

# Function to extract all up/down-regulated genes
xtract_all_deg_regTypes <- function(data = diff_exp_genes, reg_type){
  ids <- data %>%
  filter(regulation == reg_type) %>%
  select(gene_id, gene_name, log2fc) %>%
  as.data.frame() %>%
  distinct(gene_id, gene_name, .keep_all = T) %>%
  set_rownames(.$gene_id)
  
  return(ids)
}

all_deg_tables <- list(
  # Extract total background genes
  all_bg_geneIDs = crude_diff_exp_genes %>%
    as.data.frame() %>% 
    select(gene_id, gene_name, log2fc) %>%
    distinct(gene_id, gene_name, .keep_all = T) %>%
    set_rownames(.$gene_id),
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