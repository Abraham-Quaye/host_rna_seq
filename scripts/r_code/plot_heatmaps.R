#!/usr/bin/env Rscript

library(pheatmap)
library(patchwork)
library(ggplotify)

source("scripts/r_code/venn_diagram.R")
source("scripts/r_code/plot_degs.R")

## select columns to use for making heatmap
prep_heatmap_data <- function(tp){
  
  heat_data <- diff_exp_genes %>%
    filter(timepoint == tp) %>%
    arrange(qval) %>%
    slice_head(n  = 60) %>%
    arrange(desc(log2fc))
  
  r_names <- heat_data %>%
    pull(gene_name)
  
  if(tp == "t12"){
    c_names <- c("Infected 12hrs1", "Infected 12hrs2",
                 "Mock 12hrs1", "Mock 12hrs2")
  }else{
    c_names <- c("Infected 24hrs1", "Infected 24hrs2", "Infected 24hrs3",
                 "Mock 24hrs1", "Mock 24hrs2")
  }
  
  if(tp == "t12"){
    heat_data <- heat_data %>% select(-inf_fpkm_r2)
  }
  
  heat_matrix <- heat_data %>%
    select(inf_fpkm_r1:mock_fpkm_r2) %>%
    data.matrix(.) %>%
    set_rownames(r_names) %>%
    set_colnames(c_names) %>%
    t(.) %>%
    scale(.) %>%
    t(.)
  
  return(heat_matrix)
}

timpnts <- c("t12", "t24")

heatmap_data_list <- map(timpnts, prep_heatmap_data) %>%
  set_names(timpnts)

# function to plot heatmaps
plot_heatmap <- function(tp){
  tpp <- paste0(parse_number(tp), "-hpi")
  
  ht <- pheatmap(heatmap_data_list[[tp]],
                 main = paste0("Top 60 Significant DEGs at " , tpp),
                 cellwidth = 50,
                 cellheight = 10,
                 color = colorRampPalette(
                   colors = c('blue','white','red'))(250),
                 angle_col = 45,
                 fontsize = 15,
                 fontsize_row = 10,
                 fontsize_col = 15,
                 treeheight_row = 0
                 )
  if(!is.null(dev.list())){
    graphics.off()
  }
  return(ht)
}

# save heatmaps in a tibble with appropriate identifier (timepoint)
plotted_heatmaps <- tibble(timepoints = timpnts,
                           plts = map(timepoints, plot_heatmap))

p12 <- as.ggplot(plotted_heatmaps[, "plts"][[1]][[1]])
p24 <- as.ggplot(plotted_heatmaps[, "plts"][[1]][[2]])

patch_heat <- (p12 | p24) +
  plot_layout(tag_level = "new", widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "1", tag_prefix = "B") &
  theme(plot.tag = element_text(size = 22, face = "bold"))

tops <- (deg_bar_plt | patch_heat)  +
  plot_layout(widths = c(1, 1.2))

tops[[2]] <- tops[[2]] + plot_layout(tag_level = "new")
tops <- tops + plot_annotation(tag_levels = list(c("A", "B"), 1))

total_plts <- (tops / fplot) +
  plot_layout(heights = c(1.2, 1))

total_plts[[2]] <- total_plts[[2]] + plot_layout(tag_level = "new")
total_plts <- total_plts +
  plot_annotation(tag_levels = list(c("A", "B", "C"), 1)) &
  theme(plot.tag = element_text(size = 22, face = "bold"))

ggsave(plot = total_plts, filename = "results/r/figures/deg_patch_fig.png",
       width = 20, height = 20, dpi = 500)

ggsave(plot = deg_bar_plt, filename = "results/r/figures/deg_barplot.png",
       width = 8, height = 8, dpi = 400)

ggsave(plot = fplot, filename = "results/r/figures/venn_diagrams.png",
       width = 16, height = 10, dpi = 400)

ggsave(plot = p12, filename = "results/r/figures/heatmap12hpi.png",
       width = 7, height = 12, dpi = 400)

ggsave(plot = p24, filename = "results/r/figures/heatmap24hpi.png",
       width = 7, height = 12, dpi = 400)
