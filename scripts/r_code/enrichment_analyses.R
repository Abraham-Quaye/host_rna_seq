#!/usr/bin/env Rscript

library(enrichR)
library(pathview)
source("scripts/r_code/extract_deg_geneIDs.R")

dbs <- listEnrichrDbs()

active_dbs <- pull(dbs %>% filter(str_detect(libraryName, "^(Panther|GO_|KEGG)")), libraryName)

rich_down12hrs <- enrichr(pull(all_deg_tables[["up_degIDs_24hrs"]], gene_name),
                         active_dbs)

see_pathway4hrs_down <- pathview(select(all_deg_tables[["down_degIDs_12hrs"]], log2fc),
                                 pathway.id = "D08798",
                                 species = "turkey", out.suffix = "down_12hrs",
                                 )

