#!/usr/bin/env Rscript

source("scripts/r_code/extract_myDEG_geneIDs.R")

# Save geneID tables
for(t in names(all_deg_tables)){
  print(paste("Saving:", t))
  write.table(all_deg_tables[[t]]["gene_id"],
              file = paste0("results/r/tables/", t, ".txt"),
              quote = F, row.names = F, col.names = F)
}

