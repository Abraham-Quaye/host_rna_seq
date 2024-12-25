#!/usr/bin/env Rscript

library(tidyverse)

trimStats <- read.table("raw_files/ReadsQC.txt", sep = "\t", skip = 2,
                       col.names = c("sample", "raw_reads", "rawsize",
                                     "trimmed_reads", "trimmedsize", "perc_trim",
                                     "q20", "q30", "gc_content")) %>%
  select(sample, raw_reads, q20, q30, gc_content)

mapStats <- read.table("raw_files/mapped_stat.txt", header = T, sep = "\t") %>% 
  rename_with(.fn = tolower) %>%
  select(sample:multi.mapped.reads) %>%
  dplyr::rename(trimmed_reads = valid.reads) %>%
  inner_join(.,trimStats, by = "sample") %>%
  select(sample, raw_reads, everything()) %>%
  mutate(mapped_reads = str_replace(mapped.reads, "^(\\d+)\\(.*\\)$", "\\1"),
         unq_map_reads = str_replace(unique.mapped.reads, "^(\\d+)\\(.*\\)$", "\\1"),
         multi_map_reads = str_replace(multi.mapped.reads, "^(\\d+)\\(.*\\)$", "\\1")
         ) %>%
  select(sample, raw_reads, trimmed_reads, mapped_reads,
         unq_map_reads, multi_map_reads, q20, q30, gc_content) %>%
  mutate(across(raw_reads:multi_map_reads, as.numeric)) %>% 
  as_tibble() %>%
  mutate(perc_mapped = round((mapped_reads/trimmed_reads)*100, 2),
         perc_unq_mapped = round((unq_map_reads/trimmed_reads)*100, 2),
         perc_multi_mapped = round((multi_map_reads/trimmed_reads)*100, 2)) %>%
  mutate(raw_reads = round(raw_reads/1e6, 1),
         trimmed_reads = round(trimmed_reads/1e6, 1),
         mapped_reads = base::paste0(round(mapped_reads/1e6, 1), " (", perc_mapped, "%)"),
         unq_map_reads = base::paste0(round(unq_map_reads/1e6, 1), " (", perc_unq_mapped, "%)"),
         multi_map_reads = base::paste0(round(multi_map_reads/1e6, 1), " (", perc_multi_mapped, "%)")) %>%
  select(-starts_with("perc"))

# trim_rdsF <- read.table("results/fastqc_countF_Reads.txt", sep = "-", header = T)
# trim_rdsR <- read.table("results/fastqc_countR_Reads.txt", sep = "-", header = T)
# 
# trim_rdsAll <- inner_join(trim_rdsF, trim_rdsR, by = "sample") %>%
#   mutate(trimmed_reads = reads.x + reads.y,
#          sample = trimws(sample)) %>%
#   select(sample, trimmed_reads)
# 
# mapped_rdsAll <- read.csv("results/countAll_Mapped_reads.txt", header = T) %>%
#   mutate(sample = sub(pattern = "^sorted_([IU]_\\d{1,2}hrs[SN]\\d)\\.bam$", replacement = "\\1", x = bamfile)) %>%
#   select(sample, mapped_reads)
# 
# qc_stats <- tibble(sample = trimStats$sample,
#                    raw_reads = trimStats$raw_reads) %>%
#   inner_join(., trim_rdsAll, by = "sample") %>%
#   mutate(perc_ofTrimmed = round((trimmed_reads/raw_reads)*100, 1)) %>%
#   inner_join(., mapped_rdsAll, by = "sample")
# 
