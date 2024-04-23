#!/usr/bin/env Rscript

library(ggVennDiagram)
library(ggvenn)
library(patchwork)


source("scripts/r_code/extract_deg_geneIDs.R")

up <- all_deg_tables[c("up_degIDs_12hrs", "up_degIDs_24hrs")]

v_up <- list(`12hpi` = up$up_degIDs_12hrs %>% pull(gene_id),
             `24hpi` = up$up_degIDs_24hrs %>% pull(gene_id))

down <- all_deg_tables[c("down_degIDs_4hrs",
                         "down_degIDs_12hrs",
                         "down_degIDs_24hrs")]

v_down <- list(`4hpi` = down$down_degIDs_4hrs %>% pull(gene_id),
               `12hpi` = down$down_degIDs_12hrs %>% pull(gene_id),
               `24hpi` = down$down_degIDs_24hrs %>% pull(gene_id))


venn_plot <- function(input_sets){
  
  pl <- ggvenn(input_sets, stroke_size = 0.1,
               set_name_size = 8, text_size = 5,
               fill_color = c("red", "blue", "grey30", "green")) +
    coord_cartesian(clip = "off")
  
  return(pl)
}

plts <- map(list(v_up, v_down), venn_plot) %>%
  set_names(c("up", "down"))

fplot <- (plts[["up"]] |plot_spacer() | plts[["down"]]) +
  plot_layout(widths = c(1, 0.05, 1), heights = c(0.5, 0.5, 1)) +
  plot_annotation(title = "Overlap of Differentially Expressed Genes at Different Time Points",
                  tag_levels = "A",
                  theme = theme(plot.title = element_text(size = 30,
                                                          hjust = 0.5,
                                                          face = "bold",
                                                          margin = margin(b = 10)))
                  ) &
  theme(plot.tag = element_text(size = 32, face = "bold",
                                margin = margin(b = -20)))

ggsave(plot = fplot, filename = "results/r/figures/deg_vennDiagram.png",
       width = 14, height = 8, dpi = 350)
