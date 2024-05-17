#!/usr/bin/env Rscript


library(gprofiler2)
library(ggsci)

source("scripts/r_code/extract_deg_geneIDs.R")

gene_lists <- tibble(genesets = names(all_deg_tables),
                     gene_ids = map(all_deg_tables[genesets],
                                    \(df){pull(df, gene_id)}))

g_results <- gost(gene_lists$gene_ids[-1],
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
  df %>%
    ggplot(aes(reorder(term_name, recall), recall,
               group = padj, color = padj, size = gene_count)) +
    geom_point() +
    coord_flip() +
    theme_classic() +
    labs() +
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

annotate_plt <- function(plt, annot){
  
  return(plt +
    labs(title = annot[["t"]],
         x = annot[["x"]],
         y = annot[["y"]]))
}

tab_res <- tibble(timepoint = c("t4", "t12", "t24", "t72"),
                  res = list(filter(g_results, str_detect(query, "_4hrs"))
                             ,filter(g_results, str_detect(query, "_12hrs")),
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
                  plt_anot = list(list(x = "GO Term:BP",
                                       y = "Rich Factor",
                                       t = paste0("GO Enrichment: 4hpi")),
                                  list(x = "GO Term:BP",
                                       y = "Rich Factor",
                                       t = paste0("GO Enrichment: 12hpi")),
                                  list(x = "GO Term:BP",
                                       y = "Rich Factor",
                                       t = paste0("GO Enrichment: 24hpi")),
                                  list(x = "GO Term:BP",
                                       y = "Rich Factor",
                                       t = paste0("GO Enrichment: 72hpi"))),
                  p_up_bp_f = map2(p_up_bp, plt_anot, ~annotate_plt(plt = .x, annot = .y)),
                  p_down_bp_f = map2(p_down_bp, plt_anot, ~annotate_plt(plt = .x, annot = .y)),
                  plot_name = paste0("results/r/figures/go_enrich_", parse_number(timepoint))
                  )

plts <- tab_res %>% select(p_up_bp_f, p_down_bp_f, plot_name)

map2(plts$p_up_bp_f[c(2,3)], plts$plot_name[c(2,3)],
     ~ggsave(plot = .x, filename = paste0(.y, "up.png"),
             width = 14, height = 12, dpi = 400))

map2(plts$p_down_bp_f[c(2,3)], plts$plot_name[c(2,3)],
     ~ggsave(plot = .x, filename = paste0(.y, "down.png"),
             width = 14, height = 12, dpi = 400))
