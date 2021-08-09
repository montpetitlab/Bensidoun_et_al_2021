library(dplyr)
library(readr)
library(rjson)
library(purrr)
library(tibble)
setwd("~/projects/2021_bensidoun/")

fastp_multiqc <- read_tsv("outputs/fastp_trimmed/multiqc_data/multiqc_general_stats.txt") %>%
  rename_with(function(x){gsub("fastp_mqc-generalstats-fastp-", "", x)}) %>%
  select(sample = Sample, pct_duplication)

fastp_jsons<- Sys.glob("outputs/fastp_trimmed/*json") %>%
  set_names() %>%
  map_dfr(~fromJSON(file=.)$summary$before_filtering, .id = "sample") %>%
  mutate(sample = gsub("outputs\\/fastp_trimmed\\/", "", sample)) %>%
  mutate(sample = gsub("\\.trimmed\\.fastp\\.json", "", sample)) %>%
  left_join(fastp_multiqc)

counted_rate <- read_tsv("outputs/counts/raw_counts.tsv") %>%
  filter(!grepl("__", gene)) %>%
  select(!starts_with("AP_D"), -gene) %>%
  colSums() %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  select(sample, counted_reads = ".")

metadata <- left_join(fastp_jsons, counted_rate) %>%
  mutate(pct_counted = counted_reads/total_reads)


nonribo <- Sys.glob("outputs/fastp_ribo/*-nonribo*json") %>%
  set_names() %>%
  map_dfr(~fromJSON(file=.)$summary$before_filtering, .id = "sample") %>%
  mutate(sample = gsub("outputs\\/fastp_ribo\\/", "", sample)) %>%
  mutate(sample = gsub("\\-nonribo\\.fastp\\.json", "", sample)) %>%
  select(sample, non_ribo_reads = total_reads)
ribo <- Sys.glob("outputs/fastp_ribo/*-ribo*json") %>%
  set_names() %>%
  map_dfr(~fromJSON(file=.)$summary$before_filtering, .id = "sample") %>%
  mutate(sample = gsub("outputs\\/fastp_ribo\\/", "", sample)) %>%
  mutate(sample = gsub("\\-ribo\\.fastp\\.json", "", sample)) %>%
  select(sample, ribo_reads = total_reads)

metadata <- metadata %>%
  left_join(ribo) %>%
  left_join(nonribo) %>%
  mutate(pct_ribo = ribo_reads/total_reads) %>%
  mutate(pct_non_ribo = non_ribo_reads/total_reads)

flagstat <- read_tsv("outputs/star/multiqc_data/multiqc_samtools_flagstat.txt") %>%
  mutate(reads_mapped = total_passed - secondary_passed) %>%
  select(sample = Sample, reads_mapped)

metadata <- left_join(metadata, flagstat, by = "sample") %>%
  mutate(pct_mapped = reads_mapped/total_reads)

