#!/usr/bin/env Rscript

library(tidyverse)
library(readxl)
library(magrittr)
library(ggtext)
library(pheatmap)

# load in files for differentially expressed genes
diff_files <- list.files("raw_analysis/diff_exp",
                         pattern = "diff_gene_exp_\\d{1,2}hrs.xlsx",
                         full.names = T) %>%
  setNames(c("t12", "t24", "t4", "t72"))

# function to read in all the files in a compatible manner
pull_exp_data <- function(workbook){
  df <- read_excel(workbook, sheet = "Sheet1",
             col_names = T)
  
  # insert column of NA values for the third infected replicate
  if(workbook == diff_files["t12"]){
    df <- df %>%
      add_column("12" = NA, .after = 11)
  }
  
  # set the column names
  df <- df %>%
    set_colnames(c("gene_id", "gene_name", "transcript_id", "go_term",
                   "kegg", "ko_entry", "ec", "description", "trans_type",
                   "inf_fpkm_r1", "inf_fpkm_r2", "inf_fpkm_r3", "mock_fpkm_r1",
                   "mock_fpkm_r2", "fc", "log2fc", "pval", "qval", "regulation",
                   "signif"))
  
  return(df)
}

# # function to ensure consistency with groups/condition/treatment
check_grp_consistency <- function(grp_fpkms){
  # condition to reduce stringency (max(grp_fpkms, na.rm = T) > 1) & -- add to if-statement
  if(max(grp_fpkms, na.rm = T) >= (2*min(grp_fpkms, na.rm = T))){
    return(NA)
  }else{
    return(mean(grp_fpkms, na.rm = T))
  }
}

# # load all differentially expressed genes from all time points
crude_diff_exp_genes <- map_dfr(diff_files, pull_exp_data, .id = "timepoint") %>%
  mutate(inf_mean_fpkm = apply(.[, c("inf_fpkm_r1", "inf_fpkm_r2", "inf_fpkm_r3")], 1, check_grp_consistency),
         mock_mean_fpkm = apply(.[, c("mock_fpkm_r1", "mock_fpkm_r2")], 1, check_grp_consistency)) %>%
  drop_na(inf_mean_fpkm, mock_mean_fpkm) %>% 
# filter for DEGs
  filter(qval <= 0.05) %>%
  mutate(gene_known = ifelse(str_detect(gene_name, "LOC"), "unknown", "known"))

# prepare data for visualization
p_data <- crude_diff_exp_genes %>%
  group_by(timepoint, regulation) %>%
  reframe(num_genes = n()) %>%
  mutate(pos_shift = ifelse(regulation == "up", 0.2, -0.2),
         timepoint = factor(timepoint, levels = c("t4", "t12", "t24", "t72")))

# plot of differentially expressed gene tallies
deg_bar_plt <- p_data %>% 
  ggplot(aes(timepoint, num_genes, fill = regulation)) +
  geom_col(position = "dodge") +
  geom_text(aes(label = num_genes), position = position_dodge(width = 0.9),
            vjust = -0.5, fontface = "bold", size = 8) +
  scale_fill_manual(values = c("#0000ff", "#ff0000"),
                    breaks = c("down", "up"),
                    labels = c("Down Regulated", "Up Regulated")) +
  scale_y_continuous(limits = c(0, 2550),
                     expand = c(0, 0)) +
  scale_x_discrete(expand = c(0.16, 0.16),
                   breaks = c("t4", "t12", "t24"),
                   labels = c("4h.p.i", "12h.p.i", "24h.p.i")) +
  labs(title = "Differentially Expressed Genes of Turkey MDTC-RP19 Cells <br>after Infection with Turkey Hemorrhagic Enteritis Virus",
       y = "Number of Genes",
       x = element_blank(),
       fill = element_blank()) +
  theme_classic() +
  theme(plot.margin = margin(rep(20, 4)),
        plot.title = element_markdown(face = 'bold',
                                      size = 25,
                                      colour = 'black',
                                      hjust = 0.5,
                                      lineheight = -10,
                                      vjust = 0.5),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major.y = element_line(color = 'grey50',
                                          linewidth = 0.1,
                                          linetype = 1),
        panel.grid.minor = element_line(linewidth = 0.1, colour = "grey50",
                                        linetype = "dashed"),
        #adjust axis
        axis.text.x = element_markdown(size = 18, colour = 'black', face = 'bold',
                                   margin = margin(t = 10)),
        axis.text.y = element_text(size = 18, colour = 'black', face = 'bold',
                                   margin = margin(l = 10)),
        axis.title.y = element_text(size = 22,
                                    face = 'bold',
                                    color = 'black'),
        axis.ticks = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0.1, 0.8),
        legend.text = element_text(face = 'bold', 
                                   size = 18,
                                   family = 'Arial'),
        legend.margin = margin(rep(10, 4)),
        legend.key.size = unit(1.2, "cm")
  )

ggsave(plot = deg_bar_plt, filename = "results/r/deg_bar_plt.png",
       width = 15, height = 12, dpi = 400)
  

# ============================================================================
# =============================== HEATMAPS ===================================

## select columns to use for making heatmap

prep_heatmap_data <- function(tp){
  r_names <- crude_diff_exp_genes %>%
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
  
  data <- crude_diff_exp_genes %>%
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

plotted_heatmaps <- tibble(timepoints = c("t12", "t24", "t4")) %>%
  mutate(plts = map(timepoints, plot_heatmap))
      
## Save plots
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

map(pull(plotted_heatmaps, timepoints), save_heatmaps)


