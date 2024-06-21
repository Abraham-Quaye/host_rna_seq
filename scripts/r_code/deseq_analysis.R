#!/usr/bin/env Rscript

library(DESeq2)
library(apeglm)
library(ashr)
library(tidyverse)
library(magrittr)
library(rtracklayer)
library(ggplotify)
library(patchwork)
library(pheatmap)

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
         infection = factor(infection, levels = c("mock", "infected")),
         treatment = factor(treatment, levels = c(paste0("mock_", c(4, 12, 24, 72), "hrs"),
                                                  paste0("infected_", c(4, 12, 24, 72), "hrs"))
                            )
         )


## 2. Load gene matrix ====================================
gene_matrix <- read.csv("results/abundances/count_matrix/genes_count_matrix.csv",
                        header = T)

# 3. Check matching samples in experimental metadata and count matrix ==========
stopifnot(
  all(rownames(exp_metadata) %in% colnames(gene_matrix)[-1]), # checks identity/presence
  all(rownames(exp_metadata) == colnames(gene_matrix)[-1]) # checks order
  )

get_tp_matrix <- function(tp){
 gene_matrix %>%
    select(gene_id, contains(paste0("_", tp))) %>%
    set_rownames(.$gene_id) %>%
    select(-gene_id) %>%
    as.matrix(.)
}

res_to_tibble <- function(res){
    as.data.frame(res) %>%
    arrange(padj) %>%
    drop_na(padj) %>%
    mutate(gene_id = row.names(.)) %>%
    select(gene_id, everything()) %>%
    as_tibble()
}

make_norm_counts_df <- function(ncount){
    as.data.frame(ncount) %>%
    mutate(gene_id = row.names(.)) %>%
    select(gene_id, everything()) %>%
    as_tibble()
}

make_PCA <- function(dds, tp_lab){
  rld <- vst(dds, blind = F)
  pcaData <- plotPCA(rld, intgroup = "infection", returnData = T)
  pvar <- round(attr(pcaData, "percentVar") * 100)
  
  plt <- pcaData %>% 
    ggplot(aes(PC1, PC2, colour = infection)) +
    geom_point(size = 5) +
    scale_color_manual(values = c("#ff0000", "#0000ff"),
                       breaks = c("infected", "mock"),
                       labels = c("Infected", "Mock")) +
    labs(title = paste0("THEV-infected VS Uninfected PCA: ", tp_lab),
         x = paste0("PC1: ", pvar[[1]], "% Variance"),
         y = paste0("PC2: ", pvar[[2]], "% Variance")) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0),
          plot.title.position = "plot",
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 16, colour = "#000000"),
          legend.title = element_blank(),
          legend.text = element_text(size = 18, face = "bold",
                                     margin = margin(r = 30, l = 0)),
          legend.position = "top",
          legend.direction = "horizontal",
          legend.key.size = unit(1, "cm"),
          legend.margin = margin(b = -10))
  
  return(plt)
}

make_volcanoPlot <- function(res_dds, lab_tp){
  dff <- as.data.frame(res_dds) %>%
    mutate(sig = padj <= 0.05,
           reg = case_when(sig & log2FoldChange > 0 ~ "up",
                           sig & log2FoldChange < 0 ~ "down",
                           TRUE ~ "normal")) 
  
  xmin <- dff %>% drop_na(log2FoldChange, padj) %>%
    pull(log2FoldChange) %>% min(., na.rm = T)
  xmax <- dff %>% drop_na(log2FoldChange, padj) %>%
    pull(log2FoldChange) %>% max(., na.rm = T)
  
  plt <- dff %>% 
    ggplot(aes(log2FoldChange, -log10(padj), colour = reg)) +
    geom_point(alpha = 0.8, size = 2.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    scale_x_continuous(expand = c(0, 0),
                       breaks = seq(-10, 10, 1),
                       labels = seq(-10, 10, 1),
                       limits = c((xmin - 0.1), (xmax + 0.1))) +
    coord_cartesian(clip = "off") +
    labs(title = paste0("THEV-infected VS Uninfected: ", lab_tp),
         y = expression("-Log"[10]*"(P-adjusted Values)"),
         x = expression("Log"[2]*"(Fold Change)")) +
    scale_color_manual(values = c("blue", "grey", "red"),
                       breaks = c("down", "normal", "up"),
                       labels = c("Downregulated", "Not Significant", "Upregulated")) +
    theme_bw() +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 16, colour = "#000000"),
          legend.title = element_blank(),
          legend.text = element_text(size = 18, face = "bold",
                                     margin = margin(r = 30, l = 0)),
          legend.position = "inside",
          legend.position.inside = c(0.2, 0.8),
          legend.key.size = unit(1.5, "cm"),
          legend.background = element_blank(),
          legend.key = element_rect(fill = NA))
  
  return(plt)
}

# Distance Matrix
plot_sample_dists <- function(dds, lab_tp){
  rld <- rlogTransformation(dds, blind = F)
  dists <- dist(t(assay(rld)))
  dist_mat <- as.matrix(dists)
  
  if(lab_tp == "12hrs"){
    colnames(dist_mat) <- c("Infected 12hrs1", "Infected 12hrs2",
                            "Mock 12hrs1", "Mock 12hrs2")
    rownames(dist_mat) <- colnames(dist_mat)
  }else if(lab_tp == "24hrs"){
    colnames(dist_mat) <- c("Infected 24hrs1", "Infected 24hrs2", "Infected 24hrs3",
                            "Mock 24hrs1", "Mock 24hrs2")
    rownames(dist_mat) <- colnames(dist_mat)
  }else if(lab_tp == "4hrs"){
    colnames(dist_mat) <- c("Infected 4hrs1", "Infected 4hrs2", "Infected 4hrs3",
                            "Mock 4hrs1", "Mock 4hrs2")
    rownames(dist_mat) <- colnames(dist_mat)
  }else{
    colnames(dist_mat) <- c("Infected 72hrs1", "Infected 72hrs2", "Infected 72hrs3",
                            "Mock 72hrs1", "Mock 72hrs2")
    rownames(dist_mat) <- colnames(dist_mat)
  }
  
  col_heat <- colorRampPalette(hcl.colors(9, palette = "Blues"))(255)
  hmap <- pheatmap(dist_mat,
           clustering_distance_cols = dists,
           clustering_distance_rows = dists,
           color = col_heat,
           treeheight_row = 30,
           treeheight_col = 30)
  return(as.ggplot(hmap))
}

deseq_results <- tibble(timepoint = c("4hrs", "12hrs", "24hrs", "72hrs"),
              exp_data = map(timepoint, \(.x) exp_metadata %>% filter(timepoint == .x)),
              matrix = map(timepoint, get_tp_matrix),
              deseq_obj = map2(matrix, exp_data,
                               ~DESeqDataSetFromMatrix(countData = .x,
                                                       colData = .y,
                                                       design = ~infection)),
              dds_obj = map(deseq_obj, ~DESeq(.x)),
              dds_results = map(dds_obj, ~results(.x, name = "infection_infected_vs_mock",
                                                  alpha = 0.05)),
              dds_results_df = map(dds_results, res_to_tibble),
              lfc_results = map2(dds_obj, dds_results,
                                 ~lfcShrink(.x, coef = "infection_infected_vs_mock",
                                            type = "apeglm", res = .y)),
              lfc_results_tbl = map(lfc_results, res_to_tibble),
              norm_counts = map(dds_obj, ~DESeq2::counts(.x, normalized = T))) %>%
  mutate(norm_counts_df = map(norm_counts, make_norm_counts_df),
         total_res = map2(norm_counts_df, lfc_results_tbl,
                          \(.x, .y) inner_join(.x, .y, by = "gene_id") %>%
                            select(-c(baseMean, lfcSE))),
         sig_res = map(total_res, ~filter(.x, padj <= 0.05)),
         pca_plt = map2(dds_obj, timepoint, ~make_PCA(dds = .x, tp_lab = .y)),
         volcano_plt = map2(lfc_results, timepoint,
                            ~make_volcanoPlot(res_dds = .x, lab_tp = .y)),
         dist_plt = map2(dds_obj, timepoint, ~plot_sample_dists(dds = .x, lab_tp = .y))
         )

# save DEG tables
map2(deseq_results$sig_res, deseq_results$timepoint,
     ~write.csv(.x, file = paste0("results/r/tables/signif_", .y, "DEGs.csv"),
                row.names = F))

map2(deseq_results$total_res, deseq_results$timepoint,
     ~write.csv(.x, file = paste0("results/r/tables/total_", .y, "DEGs.csv"),
                row.names = F))

# Save Data QC plots
map2(deseq_results$pca_plt[c(2, 3)], deseq_results$timepoint[c(2, 3)], 
     ~ggsave(plot = .x, filename = paste0("results/r/figures/pca_", .y, ".png"),
             width = 6, height = 6, dpi = 400))

map2(deseq_results$volcano_plt[c(2, 3)], deseq_results$timepoint[c(2, 3)], 
     ~ggsave(plot = .x, filename = paste0("results/r/figures/volcano_", .y, ".png"),
             width = 8, height = 8, dpi = 400))

map2(deseq_results$dist_plt[c(2, 3)], deseq_results$timepoint[c(2, 3)], 
     ~ggsave(plot = .x, filename = paste0("results/r/figures/distPlot_", .y, ".png"),
             width = 4, height = 4, dpi = 400))

# Composite plots
patch_pca <- (deseq_results$pca_plt[[2]] + deseq_results$pca_plt[[3]]) +
  plot_layout(heights = c(1, 1), tag_level = "new")

patch_volcano <- (deseq_results$volcano_plt[[2]] | deseq_results$volcano_plt[[3]]) +
  plot_layout(widths = c(1, 1),
              tag_level = "new",
              axes = "collect") 

patch_distPlts <- (deseq_results$dist_plt[[2]] / plot_spacer()/ deseq_results$dist_plt[[3]]) +
  plot_layout(heights = c(1, 0.02, 1), tag_level = "new")


left_plts <- (patch_pca / plot_spacer() / patch_distPlts) +
  plot_layout(nrow = 3, ncol = 1, heights = c(1, 0.02, 1))

left_plts[[1]] <- left_plts[[1]] + plot_layout(tag_level = "new")
left_plts[[3]] <- left_plts[[3]] + plot_layout(tag_level = "new")

merge_plts <- (left_plts | patch_volcano) +
  plot_layout(widths = c(0.5, 1.5))

# merge_plts[[1]][[4]] <- merge_plts[[1]][[4]] + plot_layout(tag_level = "new")
merge_plts[[2]] <- merge_plts[[2]] + plot_layout(tag_level = "new")

merge_plts <- merge_plts + plot_annotation(tag_levels = c("A", "1")) &
  theme(plot.tag = element_text(size = 28, face = "bold"))

ggsave(plot = merge_plts, filename = "results/r/figures/sample_corr_figure.png",
       width = 21.5, height = 14, dpi = 400)