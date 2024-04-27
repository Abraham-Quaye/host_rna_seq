#!/usr/bin/env Rscript


library(gprofiler2)

source("scripts/r_code/extract_deg_geneIDs.R")

gene_lists <- tibble(genesets = names(all_deg_tables),
                     gene_ids = map(all_deg_tables[genesets], \(df){pull(df, gene_id)}))

g_resulsts <- gost(gene_lists$gene_ids[-1],
                   organism = "mgallopavo",
                   multi_query = F,
                   significant = T)

g_resulsts_web <- gost(gene_lists$gene_ids[-1],
                   organism = "mgallopavo",
                   multi_query = F,
                   significant = T,
                   as_short_link = T)

g_resulsts$result %>% View()
g_resulsts$result %>% filter(str_detect(query, "12hrs")) %>% View()

gostplot(g_resulsts, capped = F, interactive = F)

