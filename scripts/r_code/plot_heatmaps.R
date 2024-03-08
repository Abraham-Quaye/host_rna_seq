#!/usr/bin/env Rscript

source("scripts/r_code/deg_analysis.R")

## select columns to use for making heatmap
prep_heatmap_data <- function(tp){
  r_names <- diff_exp_genes %>%
    filter(timepoint == tp) %>%
    pull(gene_name)
  
  if(tp == "t12"){
    c_names <-  c("12hrsS1", "12hrsS2", "12hrsN1", "12hrsN2")
  }else if(tp == "t24"){
    c_names <- c("24hrsS1", "24hrsS2", "24hrsS3", "24hrsN1", "24hrsN2")
  }else if(tp == "t4"){
    c_names <- c("4hrsS1", "4hrsS2", "4hrsS3", "4hrsN1", "4hrsN2")
  }else{
    c_names <- c("72hrsS1", "72hrsS2", "72hrsS3", "72hrsN1", "72hrsN2")
  }
  
  data <- diff_exp_genes %>%
    split(.$timepoint) %>%
    set_names(c("t12", "t24", "t4"))
  
  if(tp == "t12"){
    data[[tp]] <- data[[tp]] %>% select(-inf_fpkm_r3)
  }
  
  heat <- data[[tp]] %>%
    select(inf_fpkm_r1:mock_fpkm_r2) %>%
    data.matrix(.) %>%
    set_rownames(r_names) %>%
    set_colnames(c_names) %>%
    t(.) %>%
    scale(.) %>%
    t(.)
  
  return(heat)
}

timpnts <- c("t12", "t24", "t4")

heatmap_data_list <- map(timpnts, prep_heatmap_data) %>%
  set_names(c("t12", "t24", "t4"))

# function to plot heatmaps
plot_heatmap <- function(tp){
  tpp <- paste0(parse_number(tp), "h.p.i")
  
  ht <- pheatmap(heatmap_data_list[[tp]],
                 main = paste0("DEGs of THEV-infected Turkey B-cells (" , tpp, ")"),
                 cellwidth = 60,
                 color = colorRampPalette(colors = c('blue','white','red'))(250),
                 angle_col = 45,
                 fontsize_row = 5,
                 fontsize_col = 14,
                 treeheight_row = 0,
                 # treeheight_col = 0
                 show_rownames = F
  )
  if(!is.null(dev.list())){
    graphics.off()
  }
  return(ht)
}

# save heatmaps in a tibble with appropriate identifier (timepoint)
plotted_heatmaps <- tibble(timepoints = c("t12", "t24", "t4")) %>%
  mutate(plts = map(timepoints, plot_heatmap))

# Function to save heatmaps
save_heatmaps <- function(tp){
  tpp <- paste0(str_extract(plotted_heatmaps %>%
                              filter(timepoints == tp) %>%
                              pull(timepoints),
                            "\\d+"), "hpi")
  ggsave(plot = plotted_heatmaps[plotted_heatmaps$timepoints == tp,][[2]][[1]],
         filename = paste0("results/r/deg_heatmap_", tpp, ".png"),
         width = 10, height = 15, dpi = 400)
  graphics.off()
}

## Save heatmaps
for(p in pull(plotted_heatmaps, timepoints)){
  print(paste0("Saving heatmap for timepoint ", parse_number(p), "hrs"))
  save_heatmaps(p)
}