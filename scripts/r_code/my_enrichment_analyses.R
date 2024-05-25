#!/usr/bin/env Rscript

library(gprofiler2)
library(ggsci)
library(patchwork)

source("scripts/r_code/extract_myDEG_geneIDs.R")

gene_lists <- tibble(genesets = names(all_deg_tables),
                     gene_ids = map(all_deg_tables[genesets],
                                    \(df){pull(df, gene_id)}))

g_results <- gost(gene_lists$gene_ids,
                  organism = "mgallopavo",
                  multi_query = F,
                  significant = T)


g_results <- g_results[["result"]] %>%
  filter(!str_detect(query, "^all_")) %>%
  as_tibble()

separate_feature <- function(df, feature){
  sub_df <- df %>%
    filter(source == feature) %>%
    ungroup() %>%
    select(source, term_name, gene_count = intersection_size,
           recall, padj = p_value) %>%
    arrange(padj)
  
  return(sub_df)
}

plot_goterms <- function(df){
  goTerm <- df %>% pull(source) %>% unique(.)

  df %>%
    ggplot(aes(reorder(term_name, recall), recall,
               group = padj, color = padj, size = gene_count)) +
    geom_point() +
    coord_flip() +
    theme_classic() +
    labs(x = paste0("Go Term - " , goTerm), y = "Rich Factor") +
    scale_y_continuous(breaks = c(seq(0, 1, 0.05)),
                       labels = c(seq(0, 1, 0.05))) +
    scale_color_gradient(low = "red", high = "blue", name = "Adjusted P-value") +
    scale_size(name = "Number of Genes", range = c(2, 10)) +
    theme(plot.margin = margin(rep(10, 4)),
          plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          panel.grid.major = element_line(colour = "gray", linewidth = 0.2),
          axis.title = element_text(size = 16, face = "bold", colour = "#000000"),
          axis.text = element_text(size = 10, face = "bold"),
          legend.background = element_blank(),
          legend.title = element_text(size = 12, face = "bold"),
          legend.justification = c(0, 0),
          legend.position = "inside",
          legend.position.inside = c(0.8, 0.2))
}

tab_res <- tibble(timepoint = c("t4", "t12", "t24", "t72"),
                  res = list(filter(g_results, str_detect(query, "_4hrs")),
                             filter(g_results, str_detect(query, "_12hrs")),
                             filter(g_results, str_detect(query, "_24hrs")),
                             filter(g_results, str_detect(query, "_72hrs"))),
                  down = map(res, ~ filter(.x, str_detect(query, "^down"))),
                  up = map(res, ~ filter(.x, str_detect(query, "^up"))),
                  down_mf = map(down, ~separate_feature(.x, feature = "GO:MF")),
                  down_bp = map(down, ~separate_feature(.x, feature = "GO:BP")),
                  down_cc = map(down, ~separate_feature(.x, feature = "GO:CC")),
                  down_keg = map(down, ~separate_feature(.x, feature = "KEGG")),
                  down_hp = map(down, ~separate_feature(.x, feature = "HP")),
                  up_mf = map(up, ~separate_feature(.x, feature = "GO:MF")),
                  up_bp = map(up, ~separate_feature(.x, feature = "GO:BP")),
                  up_cc = map(up, ~separate_feature(.x, feature = "GO:CC")),
                  up_keg = map(up, ~separate_feature(.x, feature = "KEGG")),
                  up_hp = map(up, ~separate_feature(.x, feature = "HP")),
                  p_up_bp = map(up_bp, plot_goterms),
                  p_down_bp = map(down_bp, plot_goterms),
                  p_up_cc = map(up_cc, plot_goterms),
                  p_down_cc = map(down_cc, plot_goterms),
                  p_up_mf = map(up_mf, plot_goterms),
                  p_down_mf = map(down_mf, plot_goterms),
                  plot_name = paste0("results/r/figures/go_enrich_", parse_number(timepoint)))

plts <- tibble(
  plts = list(tab_res$p_up_bp[[2]] + labs(title = "Upregulated GO:BP Enrichment at 12h.p.i"),
             tab_res$p_up_cc[[2]] + labs(title = "Upregulated GO:CC Enrichment at 12h.p.i"),
             tab_res$p_up_mf[[2]] + labs(title = "Upregulated GO:MF Enrichment at 12h.p.i"),
             tab_res$p_up_bp[[3]] + labs(title = "Upregulated GO:BP Enrichment at 24h.p.i"),
             tab_res$p_up_cc[[3]] + labs(title = "Upregulated GO:CC Enrichment at 24h.p.i"),
             tab_res$p_up_mf[[3]] + labs(title = "Upregulated GO:MF Enrichment at 24h.p.i"),
             tab_res$p_down_bp[[2]] + labs(title = "Downregulated GO:BP Enrichment: 12h.p.i"),
             tab_res$p_down_cc[[2]] + labs(title = "Downregulated GO:CC Enrichment: 12h.p.i"),
             tab_res$p_down_mf[[2]] + labs(title = "Downregulated GO:MF Enrichment: 12h.p.i"),
             tab_res$p_down_bp[[3]] + labs(title = "Downregulated GO:BP Enrichment: 24h.p.i"),
             tab_res$p_down_cc[[3]] + labs(title = "Downregulated GO:CC Enrichment: 24h.p.i"),
             tab_res$p_down_mf[[3]] + labs(title = "Downregulated GO:MF Enrichment: 24h.p.i")),
  plot_name = c(paste0(tab_res$plot_name[[2]], c("upBP.png", "upCC.png", "upMF.png")),
                   paste0(tab_res$plot_name[[3]], c("upBP.png", "upCC.png", "upMF.png")),
                   paste0(tab_res$plot_name[[2]], c("downBP.png", "downCC.png", "downMF.png")),
                   paste0(tab_res$plot_name[[3]], c("downBP.png", "downCC.png", "downMF.png")))
             )

# save individua plots
map2(plts$plts, plts$plot_name,
     ~ggsave(plot = .x, filename = .y,
             width = 14, height = 12, dpi = 400))

# save composite plots as one fig
total_plts <- ((plts$plts[[1]] | plts$plts[[4]]) / (plts$plts[[7]] | plts$plts[[10]])) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(plot = total_plts, filename = "results/r/figures/patch_GO_enrich.png",
       width = 30, height = 25, dpi = 350)