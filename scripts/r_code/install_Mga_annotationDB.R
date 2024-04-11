#!/usr/bin/env Rscript

library(AnnotationForge)
library(AnnotationDbi)
library(AnnotationHub)
library(biomaRt)
library(readxl)
library(tidyverse)
library(MeSHDbi)
library(meshr)

mga_info <- read_excel("raw_analysis/count_matrix/all_genes_fpkm_exp.xlsx",
                       sheet = "Sheet1", col_types = "text",
                       na = "NA") %>%
  dplyr::mutate(gene_id = str_replace(gene_id, "gene-LOC(\\d+)", "\\1")) %>%
  dplyr::select(-c(starts_with("FPKM"), starts_with("count"), trans_type)) %>%
  dplyr::rename(GID = gene_id)

mga_annotation <- list(
  gene_info = mga_info %>% dplyr::select(GID, chr:transcript_id),
  go_info = mga_info %>% dplyr::select(GID, GO) %>% tidyr::drop_na(),
  kegg_info = mga_info %>% dplyr::select(GID, KEGG) %>% tidyr::drop_na(),
  gene_description =  mga_info %>% dplyr::select(GID, Description) %>% tidyr::drop_na(),
  koEntry_info =  mga_info %>% dplyr::select(GID, KO_ENTRY) %>% tidyr::drop_na(),
  ec_info =  mga_info %>% dplyr::select(GID, EC) %>% tidyr::drop_na()
)

makeOrgPackage(gene_info = mga_annotation[["gene_info"]],
               go_info = mga_annotation[["go_info"]],
               kegg_info = mga_annotation[["kegg_info"]],
               koEntry_info = mga_annotation[["koEntry_info"]],
               ec_info = mga_annotation[["ec_info"]],
               gene_description = mga_annotation[["gene_description"]],
               version = "1.0",
               maintainer = "Abraham Quaye <quayeabraham29@gmail.com>",
               author = "Abraham Quaye <quayeabraham29@gmail.com>",
               tax_id = "9103",
               genus = "Meleagris",
               species = "gallopavo",
               outputDir = ".")

install.packages("./org.Mgallopavo.eg.db", repos = NULL)
