counts <- read_tsv("raw_counts.tsv") %>%
  filter(!grepl("__", gene)) %>%
  select(!starts_with("AP_D")) %>%
  column_to_rownames("gene")

counts5mil <- read_tsv("raw_counts_5mil.tsv") %>%
  filter(!grepl("__", gene)) %>%
  rename_with(function(x){paste0(x, "_5mil")}) %>%
  select(!starts_with("AP_D")) %>%
  column_to_rownames("gene_5mil")

all_counts <- cbind(counts, counts5mil)

cor_all_counts <- cor(all_counts)
library(corrplot)
corrplot(cor_all_counts, method="number", order = "hclust")

summary(lm(AP_A_replicat1 ~ AP_A_replicat1_5mil, data = all_counts))
