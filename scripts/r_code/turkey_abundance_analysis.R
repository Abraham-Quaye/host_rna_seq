#!/usr/bin/env Rscript

library(DESeq2)
library(apeglm)
library(tidyverse)
library(magrittr)
library(pheatmap)
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
         infection = factor(infection, levels = c("infected", "mock")),
         treatment = factor(treatment, levels = c(paste0("mock_", c(4, 12, 24, 72), "hrs"),
                                                  paste0("infected_", c(4, 12, 24, 72), "hrs"))
                            )
         )


## 2. Load gene matrix ====================================
gene_matrix <- read.csv("results/abundances/count_matrix/genes_count_matrix.csv",
                        header = T, row.names = "gene_id") %>%
  as.matrix(.)

# 3. Check matching samples in experimental metadata and count matrix ==========
stopifnot(
  all(rownames(exp_metadata) %in% colnames(gene_matrix)), # checks identity/presence
  all(rownames(exp_metadata) == colnames(gene_matrix)) # checks order
  )

# 4. Create DESeq2 object from the matrix and metadata =========================
deseq_obj <- DESeqDataSetFromMatrix(countData = gene_matrix, colData = exp_metadata,
                       design = ~ infection)

# 5. Filter low abundance genes ================================================
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we increase the speed of count modeling within DESeq2. It can also improve visualizations, as features with no information for differential expression are not plotted in dispersion plots or MA-plots.

# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 replicates per sample.
# smallest_grp_size <- 3
# rows_keep <- rowSums(counts(deseq_obj) >= 10) >= smallest_grp_size
# 
# deseq_obj <- deseq_obj[rows_keep, ]

# 6. Run differential analysis =================================================
dds <- DESeq(deseq_obj)

# 7. Explore results
resultsNames(dds)
inf_res <- results(dds, name = "infection_mock_vs_infected")
# ordered_res <- inf_res[order(inf_res$padj), ]

lfc_res <- lfcShrink(dds, coef = "infection_mock_vs_infected", type = "apeglm")

lfc_res_tab <- lfc_res %>%
  as.data.frame() %>%
  arrange(padj) %>%
  drop_na(padj) %>%
  mutate(gene_id = row.names(.)) %>%
  select(gene_id, everything()) %>%
  as_tibble()

summary(inf_res)
summary(lfc_res)

norm_counts <- counts(dds, normalized = T) %>%
  as.data.frame() %>%
  mutate(gene_id = row.names(.)) %>%
  select(gene_id, everything()) %>%
  as_tibble()

final_sig_results <- inner_join(norm_counts, lfc_res_tab, by = "gene_id") %>%
  filter(padj <= 0.05) %>% 
  arrange(padj)

# plotDispEsts(dds)
# plotMA(inf_res)
# plotMA(lfc_res)