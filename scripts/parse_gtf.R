library(rtracklayer)

gtf <- import("~/Downloads/seq_data/GCF_000146045.2_R64_genomic.gtf")


features <- gtf %>%
  as.data.frame() %>%
  filter(type == "gene")

introns <- gtf %>% 
  as.data.frame() %>%
  filter(type == "exon") %>%
  group_by(gene_id) %>% 
  filter(n() > 1) %>% 
  ungroup() %>%
  dplyr::select(gene_id)
introns <- unique(introns$gene_id)

features <- features %>%
  mutate(introns = ifelse(gene_id %in% introns, "yes", "no"))

