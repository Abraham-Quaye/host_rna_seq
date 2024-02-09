#!/usr/bin/env Rscript

library(ballgown)
library(devtools)
library(genefilter)
library(tidyverse)

# ====================================================================
# TRANSCRIPT ABUNDANCE ANALYSIS WITH DEseq2
# ====================================================================
# Hisat2 -> FeatureCounts -> DESeq2 

# make experimental data dataframe
smpl_names <- (list.dirs("results/ballgown", full.names = F) %>% .[-1])
  
exp_info <- data.frame(sample_name = smpl_names,
                       timepoint = str_extract(smpl_names, "\\d{1,2}hrs"),
                       replicate = paste0("rep", str_extract(smpl_names, "\\d$")),
                       infection = c(rep("infected", 11), rep("uninfected", 8)))


