#!/usr/bin/env Rscript

library(DESeq2)
library(apeglm)
library(tidyverse)
library(magrittr)
library(pheatmap)
library(rtracklayer)
library(readxl)
# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS WITH DEseq2
# ====================================================================
# Hisat2 -> FeatureCounts/StringTie -> DESeq2 

# 1. Make experimental data dataframe ==========================================
smpl_names <- list.dirs("results/abundances", full.names = T) %>% .[-c(1, 21)]
  
exp_metadata <- data.frame(sample_path = smpl_names,
                       timepoint = str_extract(smpl_names, "\\d{1,2}hrs"),
                       replicate = paste0("rep", str_extract(smpl_names, "\\d$")),
                       infection = c(rep("infected", 11), rep("uninfected", 8))) %>%
  mutate(sample_name = str_extract(sample_path, "abund_[IU]_\\d{1,2}hrs[SN]\\d")) %>%
  set_rownames(.$sample_name) %>%
  select(-sample_name) %>%
  mutate(timepoint = factor(timepoint, levels = c("4hrs", "12hrs", "24hrs", "72hrs")),
         infection = factor(infection, levels = c("uninfected", "infected")))

# 2. Import gene count matrix ===============================================
ncbi_gtf <- import("raw_files/annotations/turkey_genome.gtf") %>%
  as_tibble() %>%
  select(-c(seqnames:phase, gbkey, gene_biotype:transl_table)) %>%
  mutate(entrez_id = ifelse(!str_detect(gene_id, "^LOC") & !is.na(db_xref),
                          db_xref, gene_id)) %>%
  mutate(entrez_id = ifelse(str_detect(entrez_id, ":"),
                          str_replace(entrez_id, ".+:(.+)", "\\1"),
                          entrez_id)) %>%
  mutate(entrez_id = ifelse(str_detect(entrez_id, "LOC"),
                          str_replace(entrez_id, "LOC(\\d+)", "\\1"),
                          entrez_id)) %>%
  mutate(entrez_id = ifelse(str_detect(entrez_id, "^\\d+_"),
                          str_replace(entrez_id, "(\\d+)_\\d+", "\\1"),
                          entrez_id)) %>%
  select(-c(transcript_id, db_xref, gene)) %>%
  distinct(gene_id, .keep_all = T)
  
gene_matrix <- read.csv("results/abundances/count_matrix/genes_count_matrix.csv",
                        header = T) %>%
  mutate(gene_id = ifelse(str_detect(gene_id, "\\|"),
                          str_replace(gene_id, "([a-zA-Z\\d]+)\\|[a-zA-Z\\d]+", "\\1"),
                          gene_id))

gene_matrix <- left_join(gene_matrix, ncbi_gtf, by = "gene_id")

gene_matrix <- gene_matrix %>% 
  mutate(gene_id = ifelse(!is.na(entrez_id), entrez_id, gene_id)) %>%
  set_rownames(.$gene_id) %>%
  select(-c(entrez_id, gene_id)) %>%
  as.matrix()

# 3. Check matching samples in experimental metadata and count matrix ==========
stopifnot(
  all(rownames(exp_metadata) %in% colnames(gene_matrix)), # checks identity/presence
  all(rownames(exp_metadata) == colnames(gene_matrix)) # checks order
  )

# 4. Create DESeq2 object from the matrix and metadata =========================
deseq_obj <- DESeqDataSetFromMatrix(countData = gene_matrix, colData = exp_metadata,
                       design = ~ timepoint + timepoint:infection + infection)

# 5. Filter low abundance genes ================================================
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions, there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads, we reduce the memory size of the dds data object, and we increase the speed of count modeling within DESeq2. It can also improve visualizations, as features with no information for differential expression are not plotted in dispersion plots or MA-plots.

# Here we perform pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples. A recommendation for the minimal number of samples is to specify the smallest group size, e.g. here there are 3 replicates per sample.
smallest_grp_size <- 3
rows_keep <- rowSums(counts(deseq_obj) >= 10) >= smallest_grp_size

deseq_obj <- deseq_obj[rows_keep, ]

# 6. Run differential analysis =================================================
dds <- DESeq(deseq_obj)

# 7. Explore results

inf_res <- results(dds, alpha = 0.05)
# ordered_res <- inf_res[order(inf_res$padj), ]

inf_res_tab <- inf_res %>%
  as.data.frame() %>%
  arrange(padj) %>%
  drop_na(padj)

summary(inf_res)

norm_counts <- counts(dds, normalized = T) %>%
  as.data.frame()

diffgenes <- rownames(inf_res)[ which(inf_res$padj < 0.05) ]

diffcounts <- norm_counts[diffgenes, ]

# try the analysis with bowtie2 and see if geneNames are still unusable

# ==================================================================
#          Count matrix by LCSciences
# ==================================================================
lc_count_matrix <- read_excel("raw_analysis/count_matrix/gene_count_matrix.xlsx") %>%
  as.data.frame() %>%
  mutate(gene_id = ifelse(str_detect(gene_id, "^gene-LOC"),
                          str_replace(gene_id, "^gene-LOC(\\d+)$", "\\1"),
                          gene_id)) %>%
  mutate(gene_id = base::paste0("X", gene_id, "M")) %>%
  set_rownames(.$gene_id) %>%
  select(-gene_id) %>%
  mutate()
  as.matrix()

metadata_lc <- exp_metadata %>%
  mutate(sample_name = row.names(.)) %>%
  mutate(sample_name = str_replace(sample_name, "abund_([IU]_\\d{1,2}hrs[SN]\\d)", "\\1")) %>%
  set_rownames(.$sample_name) %>%
  select(-sample_name)

# 3. Check matching samples in experimental metadata and count matrix ==========
stopifnot(
  all(rownames(metadata_lc) %in% colnames(lc_count_matrix)), # checks identity/presence
  all(rownames(metadata_lc) == colnames(lc_count_matrix)) # checks order
)

# 4. Create DESeq2 object from the matrix and metadata =========================
dseq_obj <- DESeqDataSetFromMatrix(countData = lc_count_matrix, colData = metadata_lc,
                                    design = ~ timepoint + timepoint:infection + infection)


rows_keep <- rowSums(counts(dseq_obj) >= 5)

dseq_obj <- dseq_obj[rows_keep, ]

# 6. Run differential analysis =================================================
lc_dds <- DESeq(dseq_obj)

# 7. Explore results

lc_inf_res <- results(lc_dds, alpha = 0.05)
# ordered_res <- inf_res[order(inf_res$padj), ]

lc_inf_res_tab <- lc_inf_res %>%
  as.data.frame() %>%
  arrange(padj) %>%
  drop_na(padj)

summary(lc_inf_res)

norm_counts <- counts(lc_dds, normalized = T) %>%
  as.data.frame() %>%
  mutate(entrez_id = row.names(.)) %>%
  mutate(entrez_id = str_replace(entrez_id, "^X(\\d+)M\\.?\\d*$", "\\1"))





