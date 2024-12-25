#!/usr/bin/env Rscript

library(tidyverse)
library(readxl)
library(ggtext)

## Analysis Steps to obtaining gene expression levels
# 1. normalize CT-values of target genes by CT-values of reference gene by subtracting the CT-values of reference gene from all others -> d-CT
# 2. normalize the d-CT of test samples to the d-CT of calibrator sample by subtracting the d-CT of the calibrator sample from test sample d-CTs -> ddCT
# 3. calculate fold change in expression for samples with this formula:(2^-ddCT)


# Read in raw data and rename columns of interest
# Remove extra information from data -> only expression data to be analyzed
raw_expr <- read_xls("qpcr_validation/qpcr_validation.xls",
                     sheet = "Results", skip = 7, col_names = T) %>%
  select(well = Well, sample_name = `Sample Name`, gene = `Target Name`,
         ct = Cт, threshold = `Ct Threshold`, task = Task) %>%
  drop_na(well) %>%
  mutate(treatment = case_when(str_detect(sample_name, "^n") ~ "mock",
                               str_detect(sample_name, "^s") ~ "infected",
                               TRUE ~ NA_character_))

expr_data <- raw_expr %>%
  filter(gene != "tbp", gene != "actb", task != "NTC") %>%
  select(treatment, sample_name, gene, ct) %>%
  mutate(ct = as.numeric(ct))

# Extract experiment metadata (i.e, information on reference gene(s) and calibrator sample)
expr_metadata <- raw_expr %>%
  filter(is.na(gene)) %>%
  select(well, treatment, sample_name, )

ref_gene <- expr_metadata %>% filter(well == "Endogenous Control") %>%
  pull(sample_name)

ref_sample <- expr_metadata %>% filter(well == "Reference Sample") %>% pull(treatment)

# Get mean ct values for all samples and their respective genes
grp_expr <- expr_data %>%
  group_by(treatment, sample_name, gene) %>%
  summarise(mean_ct = mean(ct), 
            sd_ct = sd(ct), 
            sq_sd_ct = sd_ct^2,
            tech_reps = n(),
            .groups = "drop")

####### Group by biological replicates and plot / calculate statistics #########

# Get mean ct values for all samples and their respective genes
grp_expr <- expr_data %>%
  group_by(treatment, gene) %>%
  summarise(mean_ct = mean(ct), 
            sd_ct = sd(ct), 
            sq_sd_ct = sd_ct^2,
            tech_reps = n(),
            .groups = "drop")

# Get the mean ct values for reference gene(s)
ref_gene_ct <- grp_expr %>%
  filter(gene == ref_gene) %>%
  select(-c(gene, tech_reps)) %>%
  dplyr::rename(refG_meanCT = mean_ct, refG_sd_ct = sd_ct, refG_sq_sd_ct = sq_sd_ct)

# Calculate Delta-CT values using the reference gene mean CT values
expr_dCT_ready <- inner_join(grp_expr, ref_gene_ct, by = "treatment") %>%
  mutate(delta_ct = mean_ct - refG_meanCT,
         sd_dCT = sqrt(sq_sd_ct + refG_sq_sd_ct),
         se_dCT = sd_dCT/sqrt(tech_reps))

# Get the mean Delta-CT values for the calibrator samples
refSample_dCT <- expr_dCT_ready %>%
  filter(treatment == ref_sample) %>%
  select(gene, refS_delta_ct = delta_ct)

# Calculate Delta-delta-CT and Relative Quantities (RQs) using the calibrator sample mean dCT values
# The calculation of ΔΔCT involves subtraction of the ΔCT calibrator value. This is subtraction of an arbitrary constant, so the standard deviation of the ΔΔCT value is the same as the standard deviation of the ΔCT value. _In other words, we just use the ΔCT standard deviation values_.

expr_ddCT_ready <- inner_join(expr_dCT_ready, refSample_dCT, by = "gene") %>%
  mutate(del_delta_ct = delta_ct - refS_delta_ct,
         rq = 2^(-del_delta_ct),
         min_rq = 2^(-(del_delta_ct + se_dCT)),
         max_rq = 2^(-(del_delta_ct - se_dCT)))


# Visualize

plot_ready_data <- expr_ddCT_ready %>%
  mutate(gene = toupper(gene),
         gene = factor(gene, levels = toupper(c("apaf1", "bmf", "fadd", "madd", "pdcd4",
                                        "edem1", "ufd1", "vcp", "eif3d",
                                        "eif3m", "rpl8", "rpl10a", "gapdh"))),
         regulation = case_when(rq > 1 ~ "up",
                                rq < 1  ~ "down",
                                TRUE ~ NA_character_)) %>%
  drop_na(regulation)

plt_accessory <- tibble(x = c(0.5, 0.5, 5.4, 5.5, 5.5, 8.4, 8.5, 8.5, 12.4),
       xend = c(5.4, 0.5, 5.4, 8.4, 5.5, 8.4, 12.4, 8.5, 12.4),
       y = c(-0.5, -0.45, -0.45, -0.5, -0.45, -0.45, 0.5, 0.45, 0.45),
       yend = c(-0.5, -0.55, -0.55, -0.5, -0.55, -0.55, 0.5, 0.55, 0.55))

qpcr_plt <- plot_ready_data %>%
  ggplot(aes(gene, log2(rq), fill = regulation)) +
  geom_col() +
  geom_errorbar(aes(ymin = log2(min_rq), ymax = log2(max_rq)),
                position = position_dodge(0.9), width = 0.4) +
  geom_hline(yintercept = 0, color = "#000000") +
  annotate(geom = "text", x = c(1:12), y = c(rep(-0.1, 8), rep(0.1, 4)),
           label = levels(plot_ready_data$gene)[-13], size = 6,
           fontface = "bold.italic") +
  geom_segment(data = plt_accessory,
               aes(x = x, xend = xend, y = y, yend = yend),
               linewidth = 1.5, inherit.aes = F,
               color = "#000000") +
  annotate(geom = "text", x = c(3, 7, 10.5), y = c(rep(-0.8, 2), 0.8),
           label = paste0(c("Apoptosis\n", "ERAD\n", "Protein Synthesis\n"),
                          " Pathway"),
           size = 8, fontface = "bold") +
  scale_fill_manual(values = c("#ff0000", "#0000ff"),
                    breaks = c("up", "down"),
                    labels = c("Upregulated", "Downregulated")) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(-2.5, 3, 0.5)) +
  coord_cartesian(clip = "off", ylim = c(-2.5, NA)) +
  labs(x = element_blank(), 
       y = "Log<sub>2</sub>(Relative Quantities)",
       fill = element_blank(),
       title = "RT-qPCR Validation of Select DEGs") +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey70", linetype = "dashed",
                                          linewidth = 0.3),
        plot.title = element_markdown(size = 30, face = "bold", hjust = 0.5,
                                      colour = "#000000", margin = margin(t = 10, b = 10)),
        axis.title.y = element_markdown(size = 24, face = "bold", colour = "#000000"),
        axis.text.y = element_text(size = 16, face = "bold", colour = "#000000"),
        axis.text.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks = element_blank(),
        legend.text = element_text(size = 14, face = "bold", colour = "#000000"),
        legend.key.height = unit(0.8, "cm"),
        legend.key.spacing.x = unit(1, "cm"),
        legend.key = element_rect(fill = NA),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.90),
        legend.direction = "horizontal",
        legend.text.position = "top",
        legend.background = element_rect(colour = "gray70"))

ggsave(plot = qpcr_plt, filename = "results/r/figures/qpcr_validation.png",
       dpi = 350, width = 12, height = 10)
