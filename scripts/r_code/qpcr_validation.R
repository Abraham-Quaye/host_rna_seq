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


# Get the mean ct values for reference gene(s)
ref_gene_ct <- grp_expr %>%
  filter(gene == ref_gene) %>%
  group_by(treatment) %>%
  summarise(refG_meanCT = mean(mean_ct), refG_sd_ct = mean(sd_ct),
            refG_sq_sd_ct = mean(sq_sd_ct), .groups = "drop")


# Calculate Delta-CT values using the reference gene mean CT values
expr_dCT_ready <- inner_join(grp_expr, ref_gene_ct, by = "treatment") %>%
  mutate(delta_ct = mean_ct - refG_meanCT,
         sd_dCT = sqrt(sq_sd_ct + refG_sq_sd_ct),
         se_dCT = sd_dCT/sqrt(tech_reps))

# Get the mean Delta-CT values for the calibrator samples
refSample_dCT <- expr_dCT_ready %>%
  filter(treatment == ref_sample) %>%
  select(gene, refS_delta_ct = delta_ct) %>%
  reframe(refS_delta_ct = mean(refS_delta_ct), .by = gene)

# Calculate Delta-delta-CT and Relative Quantities (RQs) using the calibrator sample mean dCT values
# The calculation of ΔΔCT involves subtraction of the ΔCT calibrator value. This is subtraction of an arbitrary constant, so the standard deviation of the ΔΔCT value is the same as the standard deviation of the ΔCT value. _In other words, we just use the ΔCT standard deviation values_.

expr_ddCT_ready <- inner_join(expr_dCT_ready, refSample_dCT, by = "gene") %>%
  mutate(del_delta_ct = delta_ct - refS_delta_ct,
         rq = 2^(-del_delta_ct),
         min_rq = 2^(-(del_delta_ct + se_dCT)),
         max_rq = 2^(-(del_delta_ct - se_dCT))
  )

plot_ready_data <- inner_join(expr_data, expr_ddCT_ready,
                              by = c("treatment", "sample_name", "gene")) %>%
  select(treatment, sample_name, gene, se_dCT, del_delta_ct, rq, min_rq, max_rq)

# Visualize -----

all_gene_plt <- ggplot(plot_ready_data, aes(sample_name, rq, fill = gene)) +
  geom_col(position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = min_rq, ymax = max_rq),
                position = position_dodge(0.9), width = 0.4) +
  scale_fill_manual(values = rainbow(13)) +
  scale_y_continuous(expand = c(0, 0),
                     breaks = seq(0, 3, 0.25),
                     labels = format(seq(0, 3, 0.25), nsmall = 1)) +
  labs(x = "Treatment", 
       y = "Relative Quantities",
       fill = element_blank()) +
  theme_classic() +
  theme(panel.grid.major.y = element_line(colour = "grey40", linetype = "dashed",
                                          linewidth = 0.3),
        axis.title = element_text(size = 24, face = "bold", colour = "#000000"),
        axis.text.y = element_text(size = 16, face = "bold", colour = "#000000"),
        axis.text.x = element_text(size = 16, face = "bold",
                                   colour = "#000000", angle = 90),
        legend.text = element_text(size = 12, face = "bold", colour = "#000000"),
        legend.key.height = unit(0.5, "cm"),
        legend.key.spacing.x = unit(1, "cm"),
        legend.key = element_rect(fill = NA),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.text.position = "top") +
  guides(fill = guide_legend(nrow = 1))

####### Group by biological replicates and plot / calculate statistics #########

# Get mean ct values for all samples and their respective genes
grp_expr_2 <- expr_data %>%
  group_by(treatment, gene) %>%
  summarise(mean_ct = mean(ct), 
            sd_ct = sd(ct), 
            sq_sd_ct = sd_ct^2,
            tech_reps = n(),
            .groups = "drop")

# Get the mean ct values for reference gene(s)
ref_gene_ct_2 <- grp_expr_2 %>%
  filter(gene == ref_gene) %>%
  select(-c(gene, tech_reps)) %>%
  dplyr::rename(refG_meanCT = mean_ct, refG_sd_ct = sd_ct, refG_sq_sd_ct = sq_sd_ct)

# Calculate Delta-CT values using the reference gene mean CT values
expr_dCT_ready_2 <- inner_join(grp_expr_2, ref_gene_ct_2, by = "treatment") %>%
  mutate(delta_ct = mean_ct - refG_meanCT,
         sd_dCT = sqrt(sq_sd_ct + refG_sq_sd_ct),
         se_dCT = sd_dCT/sqrt(tech_reps))

# Get the mean Delta-CT values for the calibrator samples
refSample_dCT_2 <- expr_dCT_ready_2 %>%
  filter(treatment == ref_sample) %>%
  select(gene, refS_delta_ct = delta_ct)

# Calculate Delta-delta-CT and Relative Quantities (RQs) using the calibrator sample mean dCT values
# The calculation of ΔΔCT involves subtraction of the ΔCT calibrator value. This is subtraction of an arbitrary constant, so the standard deviation of the ΔΔCT value is the same as the standard deviation of the ΔCT value. _In other words, we just use the ΔCT standard deviation values_.

expr_ddCT_ready_2 <- inner_join(expr_dCT_ready_2, refSample_dCT_2, by = "gene") %>%
  mutate(del_delta_ct = delta_ct - refS_delta_ct,
         rq = 2^(-del_delta_ct),
         min_rq = 2^(-(del_delta_ct + se_dCT)),
         max_rq = 2^(-(del_delta_ct - se_dCT)))


# Visualize

plot_ready_data <- expr_ddCT_ready_2 %>%
  mutate(gene = toupper(gene),
         gene = factor(gene, levels = toupper(c("apaf1", "bmf", "fadd", "madd", "pdcd4",
                                        "edem1", "ufd1", "vcp", "eif3d",
                                        "eif3m", "rpl8", "rpl10a", "gapdh"))),
         regulation = case_when(rq > 1 ~ "up",
                                rq < 1  ~ "down",
                                TRUE ~ NA_character_)) %>%
  drop_na(regulation)
  
qpcr_plt <- plot_ready_data %>% 
  ggplot(aes(gene, log2(rq), fill = regulation)) +
  geom_col() +
  geom_errorbar(aes(ymin = log2(min_rq), ymax = log2(max_rq)),
                position = position_dodge(0.9), width = 0.4) +
  geom_hline(yintercept = 0, color = "#000000") +
  annotate(geom = "text", x = c(1:12), y = c(rep(-0.1, 8), rep(0.1, 4)),
           label = levels(plot_ready_data$gene)[-13], size = 7,
           fontface = "bold.italic") +
  scale_fill_manual(values = rainbow(2),
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
        axis.ticks.x = element_blank(),
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
  