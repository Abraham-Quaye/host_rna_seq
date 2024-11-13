#!/usr/bin/env Rscript


source("scripts/r_code/deg_analysis.R")

diff_exp_genes <- as_tibble(diff_exp_genes) %>%
  separate_wider_delim(cols = "gene_id", delim = "|",
                       names = c("gene_id", "gene_name"),
                       too_few = "align_start") %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))

# function to filter out DEG regulation types at different timepoints
xtract_geneID_tp <- function(data = diff_exp_genes, reg, tp){
  ids <- data %>%
    filter(timepoint == tp, regulation == reg) %>%
    dplyr::select(gene_id, gene_name, log2fc) %>%
    as.data.frame() %>% 
    distinct(gene_id, gene_name, .keep_all = T) %>%
    set_rownames(.$gene_id)
  
  return(ids)
}

# Function to extract all up/down-regulated genes
xtract_all_deg_regTypes <- function(data = diff_exp_genes, reg_type){
  ids <- data %>%
  filter(regulation == reg_type) %>%
  dplyr::select(gene_id, gene_name, log2fc) %>%
  as.data.frame() %>%
  distinct(gene_id, gene_name, .keep_all = T) %>%
  set_rownames(.$gene_id)
  
  return(ids)
}

all_deg_tables <- list(
  # Extract 4hr down-regulated genes
  down_degIDs_4hrs = xtract_geneID_tp(reg = "down", tp = "t4"),
  up_degIDs_4hrs = xtract_geneID_tp(reg = "up", tp = "t4"),
  # Extract 12hr up/down-regulated genes
  down_degIDs_12hrs = xtract_geneID_tp(reg = "down", tp = "t12"),
  up_degIDs_12hrs = xtract_geneID_tp(reg = "up", tp = "t12"),
  # Extract 24hr up/down-regulated genes
  down_degIDs_24hrs = xtract_geneID_tp(reg = "down", tp = "t24"),
  up_degIDs_24hrs = xtract_geneID_tp(reg = "up", tp = "t24"),
  # Extract 72hr up/down-regulated genes
  down_degIDs_72hrs = xtract_geneID_tp(reg = "down", tp = "t72"),
  # Extract all up/down-regulated genes
  all_down_degs = xtract_all_deg_regTypes(reg_type = "down"),
  all_up_degs = xtract_all_deg_regTypes(reg_type = "up")
)
