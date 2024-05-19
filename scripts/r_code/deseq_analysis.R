#!/usr/bin/env Rscript

library(DESeq2)
library(apeglm)
library(ashr)
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(readxl)

# sessionInfo()
# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS WITH DEseq2
# ====================================================================
# Hisat2 -> FeatureCounts/StringTie -> DESeq2 

# 1. Make experimental data dataframe ==========================================
smpl_names <- list.dirs("results/abundances", full.names = T) %>% .[-c(1, 21)]
  
exp_metadata <- data.frame(sample_path = smpl_names,
                       timepoint = str_extract(smpl_names, "\\d{1,2}hrs"),
                       replicate = paste0("rep", str_extract(smpl_names, "\\d$")),
                       infection = c(rep("infected", 11), rep("mock", 8))) %>%
  mutate(treatment = paste0(infection, "_", timepoint),
         sample_name = str_extract(sample_path, "abund_[IU]_\\d{1,2}hrs[SN]\\d")) %>%
  set_rownames(.$sample_name) %>%
  select(-sample_name) %>%
  mutate(timepoint = factor(timepoint, levels = c("4hrs", "12hrs", "24hrs", "72hrs")),
         infection = factor(infection, levels = c("mock", "infected")),
         treatment = factor(treatment, levels = c(paste0("mock_", c(4, 12, 24, 72), "hrs"),
                                                  paste0("infected_", c(4, 12, 24, 72), "hrs"))
                            )
         )


## 2. Load gene matrix ====================================
gene_matrix <- read.csv("results/abundances/count_matrix/genes_count_matrix.csv",
                        header = T)

# 3. Check matching samples in experimental metadata and count matrix ==========
stopifnot(
  all(rownames(exp_metadata) %in% colnames(gene_matrix)[-1]), # checks identity/presence
  all(rownames(exp_metadata) == colnames(gene_matrix)[-1]) # checks order
  )

get_tp_matrix <- function(tp){
 gene_matrix %>%
    select(gene_id, contains(paste0("_", tp))) %>%
    set_rownames(.$gene_id) %>%
    select(-gene_id) %>%
    as.matrix(.)
}

res_to_tibble <- function(res){
    as.data.frame(res) %>%
    arrange(padj) %>%
    drop_na(padj) %>%
    mutate(gene_id = row.names(.)) %>%
    select(gene_id, everything()) %>%
    as_tibble()
}

make_norm_counts_df <- function(ncount){
    as.data.frame(ncount) %>%
    mutate(gene_id = row.names(.)) %>%
    select(gene_id, everything()) %>%
    as_tibble()
}


tmp <- tibble(timepoint = c("4hrs", "12hrs", "24hrs", "72hrs"),
              exp_data = map(timepoint, \(.x) exp_metadata %>% filter(timepoint == .x)),
              matrix = map(timepoint, get_tp_matrix),
              deseq_obj = map2(matrix, exp_data,
                               ~DESeqDataSetFromMatrix(countData = .x,
                                                       colData = .y,
                                                       design = ~infection)),
              dds_obj = map(deseq_obj, ~DESeq(.x)),
              dds_results = map(dds_obj, ~results(.x, name = "infection_infected_vs_mock",
                                                  alpha = 0.05)),
              dds_results_df = map(dds_results, res_to_tibble),
              lfc_results = map2(dds_obj, dds_results,
                                 ~lfcShrink(.x, coef = "infection_infected_vs_mock",
                                            type = "ashr", res = .y)),
              lfc_results_tbl = map(lfc_results, res_to_tibble),
              norm_counts = map(dds_obj, ~DESeq2::counts(.x, normalized = T))) %>%
  mutate(norm_counts_df = map(norm_counts, make_norm_counts_df),
         total_res = map2(norm_counts_df, lfc_results_tbl,
                          ~inner_join(.x, .y, by = "gene_id")),
         sig_res = map(total_res, \(.x) filter(.x, padj <= 0.05) %>%
                         select(-c(baseMean, lfcSE))),
         maPlot_lfc = map(lfc_results, ~plotMA(.x)),
         maPlot_res = map(dds_results, ~plotMA(.x)),
         dispPlot = map(dds_obj, ~plotDispEsts(.x)))

# resultsNames(tmp$dds_obj[[1]])
map2(tmp$sig_res, tmp$timepoint,
     ~write.csv(.x, file = paste0("results/r/tables/signif_", .y, "DEGs.csv"),
                row.names = F))
