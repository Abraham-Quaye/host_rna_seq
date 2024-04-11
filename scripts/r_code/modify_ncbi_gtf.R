#!/usr/bin/env Rscript

library(tidyverse)
library(rtracklayer)


mga_ncbi_gtf <- import("raw_files/annotations/Mgallopavo_ncbi.gtf.gz") %>%
  as_tibble() %>%
  dplyr::arrange(gene_id)

sparse_geneIDs <- mga_ncbi_gtf %>%
  dplyr::filter(gene_id %in% c(paste0("ND", c(1:6, "4L")), paste0("COX", c(1:3)),
                               paste0("ATP", c(6,8)), "CYTB")) %>% 
  tidyr::drop_na(db_xref) %>%
  dplyr::mutate(db_xref = str_replace(db_xref, "^GeneID:(\\d+)$", "\\1")) %>%
  dplyr::select(gene_id, db_xref) %>%
  dplyr::arrange(gene_id)

mga_ncbi_gtf <- mga_ncbi_gtf %>% 
  dplyr::mutate(gene_id = ifelse(str_detect(gene_id, "^LOC"),
                                 str_replace(gene_id, "^LOC(\\d+)$", "\\1"),
                                 gene_id),
                gene_id = ifelse(str_detect(gene_id, "^[A-Za-z]") & !is.na(db_xref),
                                 str_replace(db_xref, "^GeneID:(\\d+)$", "\\1"),
                                 gene_id),
                gene_id = ifelse(str_detect(gene_id, "^unassigned"),
                                 product,
                                 gene_id))


for(row in 1:nrow(sparse_geneIDs)){
  mga_ncbi_gtf <- mga_ncbi_gtf %>% 
    dplyr::mutate(gene_id = case_when(gene_id == sparse_geneIDs$gene_id[row] ~
                                        sparse_geneIDs$db_xref[row],
                                      TRUE ~ gene_id))
}

export(mga_ncbi_gtf, "raw_files/annotations/mod_turkey_genome.gtf")


# mga_go_ncbi <- read_tsv("raw_files/annotations/Mgallopavo_GOncbi.gaf.gz",
#                         skip = 8) %>%
#   dplyr::select(c(gene_id = GeneID, symbol = Symbol, full_name = Gene_Name,
#                   go_id = GO_ID, go_evidence_code = Evidence_Code,
#                   ncbi_ref = Reference)) %>%
#   dplyr::mutate(gene_id = as.character(gene_id)) %>%
#   dplyr::arrange(gene_id)