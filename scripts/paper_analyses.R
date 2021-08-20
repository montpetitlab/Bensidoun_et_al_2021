library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(ggplot2)
library(DESeq2)
library(UpSetR)
library(xlsx)
library(rtracklayer)
library(ggpmisc)
library(ggpubr)
library(broom)

#library(clusterProfiler)

setwd("~/github/2021-bensidoun/")

# functions ---------------------------------------------------------------

orf_to_common <- function(keys){
  commons <- AnnotationDbi::select(org.Sc.sgd.db::org.Sc.sgd.db,
                                   keys = keys,
                                   columns=c("COMMON"),
                                   keytype="ORF")
  commons <- commons[match(unique(commons$ORF), commons$ORF), ]
  commons <- commons[ , -2]
  commons$COMMON <- ifelse(is.na(commons$COMMON), commons$ORF, commons$COMMON)
  commons <- unique(commons[ , ])
  commons <- commons$COMMON
  return(commons)
}

orf_to_description <- function(keys){
  description <- AnnotationDbi::select(org.Sc.sgd.db::org.Sc.sgd.db,
                                       keys = keys,
                                       columns=c("DESCRIPTION"),
                                       keytype="ORF")
  description <- description[match(unique(description$ORF), description$ORF), ]
  description <- description$DESCRIPTION
  return(description)
}

# differential expression against the total transcriptome control --------------
info <- read_csv("inputs/info.csv")
info$condition<- factor(info$condition, levels = c("input", "basket", "basketless", 
                                                   "total", "igg", "gbp"))
counts <- read_tsv("outputs/counts/raw_counts.tsv") %>%
  filter(!grepl("__", gene)) %>%
  select(!starts_with("AP_D")) %>%
  column_to_rownames("gene")

dds <- DESeqDataSetFromMatrix(counts,
                              colData = info,
                              design = ~ condition)
ds <- DESeq(dds, test="Wald")
#resultsNames(ds)
# log2 fold change (MLE): condition total vs input
# results are log2(total/input)
# log2FC values are expressed as that change in gene expression in total compared
# to that in input.

res_basket <- results(ds, contrast = c("condition", "basket", "input"), alpha = .05)
res_basket <- res_basket %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 1) %>%
  mutate(common = orf_to_common(gene)) %>%
  mutate(description = orf_to_description(gene))

res_basketless <- results(ds, contrast = c("condition", "basketless", "input"), alpha = .05)
res_basketless <- res_basketless %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 1) %>%
  mutate(common = orf_to_common(gene)) %>%
  mutate(description = orf_to_description(gene))

res_total <- results(ds, contrast = c("condition", "total", "input"), alpha = .05)
res_total <- res_total %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 1)  %>%
  mutate(common = orf_to_common(gene)) %>%
  mutate(description = orf_to_description(gene))

res_gbp <- results(ds, contrast = c("condition", "gbp", "input"), alpha = .05)
res_gbp <- res_gbp %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 1)  %>%
  mutate(common = orf_to_common(gene)) %>%
  mutate(description = orf_to_description(gene))

res_igg <- results(ds, contrast = c("condition", "igg", "input"), alpha = .05)
res_igg <- res_igg %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  filter(padj < 0.05) %>%
  filter(abs(log2FoldChange) > 1)  %>%
  mutate(common = orf_to_common(gene)) %>%
  mutate(description = orf_to_description(gene))

# write.xlsx(res_total, file="outputs/xlsx/deseq2_results_background.xlsx", sheetName="res_total", row.names=FALSE)
# write.xlsx(res_basket, file="outputs/xlsx/deseq2_results_background.xlsx", sheetName="res_basket", row.names=FALSE, append = T)
# write.xlsx(res_basketless, file="outputs/xlsx/deseq2_results_background.xlsx", sheetName="res_basketless", row.names=FALSE, append = T)
# write.xlsx(res_igg, file="outputs/xlsx/deseq2_results_background.xlsx", sheetName="res_igg", row.names=FALSE, append = T)
# write.xlsx(res_gbp, file="outputs/xlsx/deseq2_results_background.xlsx", sheetName="res_gbp", row.names=FALSE, append = T)

# remove specific background  --------------------------------------------------

remove_igg_background <- function(res, res_igg){
  res_igg_up <- res_igg %>%
    filter(log2FoldChange > 1) %>%
    filter(padj < 0.05)
  
  res_igg_down <- res_igg %>%
    filter(log2FoldChange < -1) %>%
    filter(padj < 0.05)
  
  background_up <- c(res_igg_up$gene)
  background_down <- c(res_igg_down$gene)
  
  res_up <- res %>%
    filter(log2FoldChange > 1) %>%
    filter(!gene %in% background_up)
  
  res_down <- res %>%
    filter(log2FoldChange < -1) %>%
    filter(!gene %in% background_down)
  
  res_no_background <- rbind(res_up, res_down)
  return(res_no_background)
}

remove_gbp_background <- function(res, res_gbp){
  res_gbp_up <- res_gbp %>%
    filter(log2FoldChange > 1) %>%
    filter(padj < 0.05)
  
  res_gbp_down <- res_gbp %>%
    filter(log2FoldChange < -1) %>%
    filter(padj < 0.05)
  
  background_up <- c(res_gbp_up$gene)
  background_down <- c(res_gbp_down$gene)
  
  res_up <- res %>%
    filter(log2FoldChange > 1) %>%
    filter(!gene %in% background_up)
  
  res_down <- res %>%
    filter(log2FoldChange < -1) %>%
    filter(!gene %in% background_down)
  
  res_no_background <- rbind(res_up, res_down)
  return(res_no_background)
}

res_basket_nigg <- remove_igg_background(res_basket, res_igg)
res_total_ngbp <- remove_gbp_background(res_total, res_gbp)
res_basketless_ngbp <- remove_gbp_background(res_basketless, res_gbp)

# write.xlsx(res_total_ngbp, file="outputs/xlsx/deseq2_results_no_specific_background.xlsx", sheetName="res_total_ngbp", row.names=FALSE)
# write.xlsx(res_basket_nigg, file="outputs/xlsx/deseq2_results_no_specific_background.xlsx", sheetName="res_basket_nigg", row.names=FALSE, append = T)
# write.xlsx(res_basketless_ngbp, file="outputs/xlsx/deseq2_results_no_specific_background.xlsx", sheetName="res_basketless_ngbp", row.names=FALSE, append = T)

# calculate enriched/etc. set size
res_basket_nigg %>%
  mutate(enriched = ifelse(log2FoldChange > 1, "enriched", "not_enriched")) %>%
  group_by(enriched) %>%
  tally()

res_basketless_ngbp %>%
  mutate(enriched = ifelse(log2FoldChange > 1, "enriched", "not_enriched")) %>%
  group_by(enriched) %>%
  tally()

res_total_ngbp %>%
  mutate(enriched = ifelse(log2FoldChange > 1, "enriched", "not_enriched")) %>%
  group_by(enriched) %>%
  tally()

# upset plot after removal (nsb - no specific background) ----------------------
up_nsb <- fromList(list("total" = filter(res_total_ngbp, log2FoldChange > 1)$gene,
                        'basket' = filter(res_basket_nigg, log2FoldChange > 1)$gene,
                        'basketless' = filter(res_basketless_ngbp, log2FoldChange > 1)$gene))
pdf("tmp_up_nsb_upset.pdf", width = 5, 2.5)
upset(up_nsb, nsets = 3, set_size.show = T, set_size.scale_max = 1800)
dev.off()

down_nsb <-  fromList(list("total" = filter(res_total_ngbp, log2FoldChange < -1)$gene,
                           'basket' = filter(res_basket_nigg, log2FoldChange < -1)$gene,
                           'basketless' = filter(res_basketless_ngbp, log2FoldChange < -1)$gene))
pdf("tmp_down_nsb_upset.pdf", width = 5, 2.5)
upset(down_nsb, nsets = 3, set_size.show = T, set_size.scale_max = 1800)
dev.off()


# separate genes from 7 permutations, no specific background -------------------
fromList2 <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

get_intersect_members <- function (x, ...){
  require(dplyr)
  require(tibble)
  x <- x[,sapply(x, is.numeric)][,0<=colMeans(x[,sapply(x, is.numeric)],na.rm=T) & colMeans(x[,sapply(x, is.numeric)],na.rm=T)<=1]
  n <- names(x)
  x %>% rownames_to_column() -> x
  l <- c(...)
  a <- intersect(names(x), l)
  ar <- vector('list',length(n)+1)
  ar[[1]] <- x
  i=2
  for (item in n) {
    if (item %in% a){
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '>= 1')
        i <- i + 1
      }
    } else {
      if (class(x[[item]])=='integer'){
        ar[[i]] <- paste(item, '== 0')
        i <- i + 1
      }
    }
  }
  do.call(filter_, ar) %>% column_to_rownames() -> x
  return(x)
}

format_intersect_members <- function(intersect_members_df){
  intersect_members_df <- intersect_members_df %>%
    rownames_to_column("gene") %>%
    mutate(common = orf_to_common(gene)) %>%
    mutate(description = orf_to_description(gene))
  return(intersect_members_df)
}

# res total minus background
# res basket minus background
# res basketless minus background
# res basket, basketless, total minus background
# res basket, total, minus background
# res basketless, total, minus background
# res basket, basketless, minus background

up_nsb2 <- fromList2(list("total" = filter(res_total_ngbp, log2FoldChange > 1)$gene,
                          'basket' = filter(res_basket_nigg, log2FoldChange > 1)$gene,
                          'basketless' = filter(res_basketless_ngbp, log2FoldChange > 1)$gene))

total_basket_nsb_up <- format_intersect_members(get_intersect_members(up_nsb2, "total", "basket"))
total_basketless_nsb_up <- format_intersect_members(get_intersect_members(up_nsb2, "total", "basketless"))
total_basketless_basket_nsb_up <- format_intersect_members(get_intersect_members(up_nsb2, "total", "basket", "basketless"))
basket_nsb_up <- format_intersect_members(get_intersect_members(up_nsb2, "basket"))
basketless_nsb_up <- format_intersect_members(get_intersect_members(up_nsb2, "basketless"))
total_nsb_up <- format_intersect_members(get_intersect_members(up_nsb2, "total"))
basket_basketless_nsb_up <- format_intersect_members(get_intersect_members(up_nsb2, "basket", "basketless"))

# write.xlsx(total_basket_nsb_up, file="outputs/xlsx/upset_permutations_up_nsb.xlsx", sheetName="total_basket_nsb_up", row.names=FALSE)
# write.xlsx(total_basketless_nsb_up, file="outputs/xlsx/upset_permutations_up_nsb.xlsx", sheetName="total_basketless_nsb_up", append = T, row.names=FALSE)
# write.xlsx(total_basketless_basket_nsb_up, file="outputs/xlsx/upset_permutations_up_nsb.xlsx", sheetName="total_basketless_basket_nsb_up", append = T, row.names=FALSE)
# write.xlsx(basket_nsb_up, file="outputs/xlsx/upset_permutations_up_nsb.xlsx", sheetName="basket_nsb_up", append = T, row.names=FALSE)
# write.xlsx(basketless_nsb_up, file="outputs/xlsx/upset_permutations_up_nsb.xlsx", sheetName="basketless_nsb_up", append = T, row.names=FALSE)
# write.xlsx(total_nsb_up, file="outputs/xlsx/upset_permutations_up_nsb.xlsx", sheetName="total_nsb_up", append = T, row.names=FALSE)
# write.xlsx(basket_basketless_nsb_up, file="outputs/xlsx/upset_permutations_up_nsb.xlsx", sheetName="basket_basketless_nsb_up", append = T, row.names=FALSE)

down_nsb2 <- fromList2(list("total" = filter(res_total_ngbp, log2FoldChange < -1)$gene,
                            'basket' = filter(res_basket_nigg, log2FoldChange < -1)$gene,
                            'basketless' = filter(res_basketless_ngbp, log2FoldChange < -1)$gene))

total_basket_nsb_down <- format_intersect_members(get_intersect_members(down_nsb2, "total", "basket"))
total_basketless_nsb_down <- format_intersect_members(get_intersect_members(down_nsb2, "total", "basketless"))
total_basketless_basket_nsb_down <- format_intersect_members(get_intersect_members(down_nsb2, "total", "basket", "basketless"))
basket_nsb_down <- format_intersect_members(get_intersect_members(down_nsb2, "basket"))
basketless_nsb_down <- format_intersect_members(get_intersect_members(down_nsb2, "basketless"))
total_nsb_down <- format_intersect_members(get_intersect_members(down_nsb2, "total"))
basket_basketless_nsb_down <- format_intersect_members(get_intersect_members(down_nsb2, "basket", "basketless"))

# write.xlsx(total_basket_nsb_down, file="outputs/xlsx/upset_permutations_down_nsb.xlsx", sheetName="total_basket_nsb_down", row.names=FALSE)
# write.xlsx(total_basketless_nsb_down, file="outputs/xlsx/upset_permutations_down_nsb.xlsx", sheetName="total_basketless_nsb_down", append = T, row.names=FALSE)
# write.xlsx(total_basketless_basket_nsb_down, file="outputs/xlsx/upset_permutations_down_nsb.xlsx", sheetName="total_basketless_basket_nsb_down", append = T, row.names=FALSE)
# write.xlsx(basket_nsb_down, file="outputs/xlsx/upset_permutations_down_nsb.xlsx", sheetName="basket_nsb_down", append = T, row.names=FALSE)
# write.xlsx(basketless_nsb_down, file="outputs/xlsx/upset_permutations_down_nsb.xlsx", sheetName="basketless_nsb_down", append = T, row.names=FALSE)
# write.xlsx(total_nsb_down, file="outputs/xlsx/upset_permutations_down_nsb.xlsx", sheetName="total_nsb_down", append = T, row.names=FALSE)
# write.xlsx(basket_basketless_nsb_down, file="outputs/xlsx/upset_permutations_down_nsb.xlsx", sheetName="basket_basketless_nsb_down", append = T, row.names=FALSE)


# correlations with gene features -----------------------------------------

gtf <- import("GCF_000146045.2_R64_genomic.gtf")

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

features <- features %>%
  dplyr::select(gene_id, width, strand, pseudo, introns)

norm_dds <- estimateSizeFactors(dds)
norm_counts <- DESeq2::counts(norm_dds, normalized=T)
polya_counts <- norm_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene, starts_with("AP_G")) %>%
  rowwise() %>% 
  mutate(mean_polya = mean(c(AP_G_replicat1, AP_G_replicat2, AP_G_replicat3))) %>%
  select(gene, mean_polya)

vsd <- vst(dds)
vsd <- assay(vsd)
polya_vsd <- vsd %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  select(gene, starts_with("AP_G")) %>%
  rowwise() %>% 
  mutate(mean_polya_vst = mean(c(AP_G_replicat1, AP_G_replicat2, AP_G_replicat3))) %>%
  select(gene, mean_polya_vst)

features <- features %>%
  left_join(polya_counts, by = c("gene_id" = "gene")) %>%
  left_join(polya_vsd, by = c("gene_id" = "gene"))

# format for and perform stats
# lm() -- probably not the best test to do
# wilcox.test() -- wilcoxon rank sum test. Is mean gene length different
#                  between two groups? (data are not normally distributed)
# chisq.test() -- chi squared test. Are genes with introns or from one strand
#                 more or less likely to be associated with NPCs?

## IGG
res_igg_features <- res_igg %>%
  left_join(features, by = c("gene" = "gene_id")) %>%
  mutate(change = ifelse(log2FoldChange > 0 , "log2fc_pos", "log2fc_neg"))

#summary(lm(log2FoldChange ~ width, data = res_igg_features))

pairwise.wilcox.test(res_igg_features$width, res_igg_features$change) %>%
  broom::tidy() 
# pairwise.wilcox.test: Is mean gene length different between two groups? 
#   group1     group2        p.value
#  log2fc_pos log2fc_neg 0.00000124
# Yes, mean gene length is different between genes that are enriched in 
# res_igg over controls
chisq.test(table(res_igg_features$change, res_igg_features$introns))
# genes with introns are not more likely to be enriched in igg
chisq.test(table(res_igg_features$change, as.character(res_igg_features$strand)))
# genes from a one strand are not more likely to be enriched in igg

## GBP
res_gbp_features <- res_gbp %>%
  left_join(features, by = c("gene" = "gene_id")) %>%
  mutate(change = ifelse(log2FoldChange > 0 , "log2fc_pos", "log2fc_neg"))

# summary(lm(log2FoldChange ~ width, data = res_gbp_features))
pairwise.wilcox.test(res_gbp_features$width, res_gbp_features$change) %>%
  broom::tidy()
# pairwise.wilcox.test: Is mean gene length different between two groups? 
#   group1     group2        p.value
#   log2fc_pos log2fc_neg    5.23e-14
# Yes, mean gene length is different between genes that are enriched in 
# res_gbp over controls
chisq.test(table(res_gbp_features$change, res_gbp_features$introns))
# genes with introns are not more likely to be enriched in gbp
chisq.test(table(res_gbp_features$change, as.character(res_gbp_features$strand)))
# genes from a one strand are not more likely to be enriched in gbp

## BASKET with specific backgroun (IGG) subtracted 
res_basket_features <- res_basket_nigg %>%
  left_join(features, by = c("gene" = "gene_id")) %>%
  mutate(change = ifelse(log2FoldChange > 0 , "log2fc_pos", "log2fc_neg"))

#summary(lm(log2FoldChange ~ width, data = res_basket_features))
pairwise.wilcox.test(res_basket_features$width, res_basket_features$change) %>%
  broom::tidy()
# Yes, mean gene length is different between genes that are enriched in 
# res_basket_nigg over controls
#   group1     group2       p.value
#   log2fc_pos log2fc_neg   5.71e-235
chisq.test(table(res_basket_features$change, res_basket_features$introns))
# genes with introns are less likely to be enriched in basket (p-value < 2.2e-16)
chisq.test(table(res_basket_features$change, as.character(res_basket_features$strand)))
# genes from a one strand are not more likely to be enriched in basket

## BASKETLESS with specific background removed (GBP)
res_basketless_features <- res_basketless_ngbp %>%
  left_join(features, by = c("gene" = "gene_id")) %>%
  mutate(change = ifelse(log2FoldChange > 0 , "log2fc_pos", "log2fc_neg"))

#summary(lm(log2FoldChange ~ width, data = res_basketless_features))
pairwise.wilcox.test(res_basketless_features$width, res_basketless_features$change) %>%
  broom::tidy()
# Yes, mean gene length is different between genes that are enriched in 
# res_basket_nigg over controls
#   group1     group2       p.value
#   log2fc_pos log2fc_neg   3.36e-125
chisq.test(table(res_basketless_features$change, res_basketless_features$introns))
# genes with introns are less likely to be enriched in basketless (p-value < 2.2e-16)
chisq.test(table(res_basketless_features$change, as.character(res_basketless_features$strand)))
# genes from a one strand are not more likely to be enriched in basketless

## TOTAL with specific background removed (GBP)
res_total_features <- res_total_ngbp %>%
  left_join(features, by = c("gene" = "gene_id")) %>%
  mutate(change = ifelse(log2FoldChange > 0 , "log2fc_pos", "log2fc_neg"))

#summary(lm(log2FoldChange ~ width, data = res_total_features))
pairwise.wilcox.test(res_total_features$width, res_total_features$change) %>%
  broom::tidy()
# Yes, mean gene length is different between genes that are enriched in 
# res_basket_nigg over controls
#   group1     group2       p.value
#   log2fc_pos log2fc_neg   1.39e-173
chisq.test(table(res_total_features$change, res_total_features$introns))
# genes with introns are less likely to be enriched in total (p-value < 2.2e-16)
chisq.test(table(res_total_features$change, as.character(res_total_features$strand)))
# genes from a one strand are not more likely to be enriched in total

# plots -------------------------------------------------------------------
lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

my_formula <- y ~ x


plt1 <- ggplot(res_basketless_features, aes(x = width)) +
  geom_density(alpha = .5, aes(fill = change)) +
  geom_density(data = features, aes(x = width), linetype = "dashed", color ="grey") +
  theme_minimal() +
  ggtitle("Basketless - GFP")+ 
  xlim(0, 15000) +
  ylim(0, 0.0015)

plt2 <- ggplot(res_total_features, aes(x = width)) +
  geom_density(alpha = .5, aes(fill = change)) +
  geom_density(data = features, aes(x = width), linetype = "dashed", color ="grey") +
  theme_minimal() +
  ggtitle("Total - GFP")+ 
  xlim(0, 15000) +
  ylim(0, 0.0015)

plt3 <- ggplot(res_basket_features, aes(x = width)) +
  geom_density(alpha = .5, aes(fill = change)) +
  geom_density(data = features, aes(x = width), linetype = "dashed", color ="grey") +
  theme_minimal() +
  ggtitle("Basket - IgG")+ 
  xlim(0, 15000) +
  ylim(0, 0.0015)

plt4 <- ggplot(res_igg_features, aes(x = width)) +
  geom_density(alpha = .5, aes(fill = change)) +
  geom_density(data = features, aes(x = width), linetype = "dashed", color ="grey") +
  theme_minimal() +
  ggtitle("IgG")+ 
  xlim(0, 15000) +
  ylim(0, 0.0015)

plt5 <- ggplot(res_gbp_features, aes(x = width)) +
  geom_density(alpha = .5, aes(fill = change)) +
  geom_density(data = features, aes(x = width), linetype = "dashed", color ="grey") +
  theme_minimal() +
  ggtitle("GFP") + 
  xlim(0, 15000) +
  ylim(0, 0.0015)

ggarrange1 <- ggarrange(plt1, plt2, plt3, nrow = 3, legend = "none", common.legend = T)
ggarrange2 <- ggarrange(plt4, plt5, common.legend = T, nrow =2,
                        legend = "bottom")

ggarrange(ggarrange1, ggarrange2, ncol = 2, common.legend = T)


# chi square plts (introns) ---------------------------------------------------

plt1 <- ggplot(res_total_features) +
  aes(x = change, fill = introns) +
  geom_bar(stat = "count") +
  geom_text(aes(label = ..count..), stat = "count") +
  scale_fill_hue() +
  theme_minimal() +
  ggtitle("Total - GFP")

plt2 <- ggplot(res_basketless_features) +
  aes(x = change, fill = introns) +
  geom_bar(stat = "count") +
  geom_text(aes(label = ..count..), stat = "count") +
  scale_fill_hue() +
  theme_minimal() +
  ggtitle("Basketless - GFP")

plt3 <- ggplot(res_basket_features) +
  aes(x = change, fill = introns) +
  geom_bar(stat = "count") +
  geom_text(aes(label = ..count..), stat = "count") +
  scale_fill_hue() +
  theme_minimal() +
  ggtitle("Basket - IgG")


plt4 <- ggplot(res_igg_features) +
  aes(x = change, fill = introns) +
  geom_bar(stat = "count") +
  geom_text(aes(label = ..count..), stat = "count") +
  scale_fill_hue() +
  theme_minimal() +
  ggtitle("IgG")

plt5 <- ggplot(res_gbp_features) +
  aes(x = change, fill = introns) +
  geom_bar(stat = "count") +
  geom_text(aes(label = ..count..), stat = "count") +
  scale_fill_hue() +
  theme_minimal() +
  ggtitle("GFP")

ggarrange1 <- ggarrange(plt1, plt2, plt3, nrow = 3, legend = "none", common.legend = T)
ggarrange2 <- ggarrange(plt4, plt5, common.legend = T, nrow =2,
                        legend = "bottom")

ggarrange(ggarrange1, ggarrange2, ncol = 2, common.legend = T)


# exploring some more. -----------------------------------------------

# FOR THINGS THAT ARE SIGNIFCANTLY ASSOCIATED WITH BASKET AND BASKETLESS,
# ARE THEY MORE LIKELY TO BE ASSOCIATED WITH BASKET THAN BASKETLESS? -----------
# Yes, on average, the l2fc is 1.2 lower for basketless

res_basket_nigg_shared <- res_basket_nigg %>%
  filter(gene %in% c(total_basketless_basket_nsb_up$gene, basket_basketless_nsb_up$gene)) %>%
  mutate(ap = "basket") 

res_basketless_ngbp_shared <- res_basketless_ngbp %>%
  filter(gene %in% c(total_basketless_basket_nsb_up$gene, basket_basketless_nsb_up$gene)) %>%
  mutate(ap = "basketless")

shared_up <- left_join(res_basket_nigg_shared, res_basketless_ngbp_shared, by = c("gene", "baseMean", "common", "description"))

ggplot(shared_up, aes(x = log2FoldChange.x, y = log2FoldChange.y, label = common)) +
  geom_point(alpha = .1) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(plot.title.position =  "plot",
        plot.title =element_text(size = 9)) +
  ylim(1, 5) +
  labs(x = "log2FC basket", y = "log2FC basketless", title = "Basket + Basketless (n = 3 + 436)")

t.test(x = shared$log2FoldChange.x, y = shared$log2FoldChange.y, paired = T)


res_basket_nigg_shared <- res_basket_nigg %>%
  filter(gene %in% c(total_basketless_basket_nsb_down$gene, basket_basketless_nsb_down$gene)) %>%
  mutate(ap = "basket") 

res_basketless_ngbp_shared <- res_basketless_ngbp %>%
  filter(gene %in% c(total_basketless_basket_nsb_down$gene, basket_basketless_nsb_down$gene)) %>%
  mutate(ap = "basketless")

shared_down <- left_join(res_basket_nigg_shared, res_basketless_ngbp_shared, by = c("gene", "baseMean", "common", "description"))

ggplot(shared_down, aes(x = log2FoldChange.x, y = log2FoldChange.y, label = common)) +
  geom_point(alpha = .1) +
  geom_abline(slope = 1, intercept = 0) +
  theme_minimal() +
  theme(plot.title.position =  "plot",
        plot.title =element_text(size = 9)) +
  ylim(-6, -1) +
  labs(x = "log2FC basket", y = "log2FC basketless", title = "Basket + Basketless (n = 43 + 459)")

# WHAT DIFFERENTIATES THESE GENES? WHAT FEATURES MAKE A GENE MORE LIKELY TO END 
# UP IN A PORE, OR IN A BASKET OR BASKETLESS PORE? -----------------------

# sets -- anything that's more likely to be in basketless vs anything thanks only
# in basket.
# basketless = 5 + 40 + 3 + 436
# basket = 746


# EXPRESSION ------

# it looks like transcripts that are more likely to be in basketless are 
# intron-less and highly expressed
# The opposite is true for transcripts that are only significantly enriched in basket 

ggplot(res_basketless_features, aes(x = log2FoldChange, y = mean_polya_vst, color = introns)) +
  geom_point(alpha = 0.25)


# stats
plt1_exp <- ggplot(res_basketless_features %>%
         filter(introns == "no"), 
       aes(x = mean_polya_vst, fill = change)) +
  geom_density(alpha = .5) +
  theme_minimal()+
  theme_minimal() +
  theme(plot.title = element_text(size = 9),
        plot.title.position = "plot") +
  labs(title = "Basketless (n = 5 + 40 + 3 + 436 +\n30 + 71 + 43 + 459)",
       x = "mean VST count in poly(A) library")

res_basketless_features %>%
  filter(introns == "no") %>%
  group_by(change) %>%
  summarize(mean_mean_polya_vst = mean(mean_polya_vst))

res_basketless_features %>%
  group_by(introns, change) %>% 
  summarise(mean_polya_vst = mean(mean_polya_vst),
            mean_polya_counts = mean(mean_polya))

res_basketless_features %>% 
  group_by(introns) %>% 
  do(tidy(t.test(mean_polya_vst ~ change, data = .)))

## Look at distinctly enriched in basket
ggplot(res_basket_features %>%
         filter(gene %in% c(basket_nsb_up$gene, basket_nsb_down$gene)), 
       aes(x = log2FoldChange, y = mean_polya_vst, color = introns)) +
  geom_point(alpha = 0.25) 

plt2_exp <- ggplot(res_basket_features %>%
         filter(gene %in% c(basket_nsb_up$gene, basket_nsb_down$gene)) %>%
         filter(introns == "no"), 
       aes(x = mean_polya_vst, fill = change)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  theme(plot.title = element_text(size = 9),
        plot.title.position = "plot") +
  labs(title = "Basket (n = 746 + 743)", 
       x = "mean VST count in poly(A) library")

res_basket_features %>%
  filter(gene %in% c(basket_nsb_up$gene, basket_nsb_down$gene)) %>%
  group_by(introns) %>% 
  do(tidy(t.test(mean_polya_vst ~ change, data = .)))

res_basket_features %>%
  filter(gene %in% c(basket_nsb_up$gene, basket_nsb_down$gene)) %>%
  group_by(introns, change) %>% 
  summarise(mean_polya_vst = mean(mean_polya_vst),
            mean_polya_counts = mean(mean_polya))

ggarrange(plt2_exp, plt1_exp, common.legend = T, nrow = 1, legend = "bottom")

# LENGTH ---------

plt1 <- ggplot(res_basketless_features %>%
                filter(gene %in% c(total_basketless_basket_nsb_up$gene, 
                                   total_basketless_basket_nsb_down$gene)), 
              aes(x = width)) +
  geom_density(alpha = .5, aes(fill = change)) +
  geom_density(data = features, aes(x = width), linetype = "dashed", color ="grey") +
  theme_minimal() +
  theme(plot.title = element_text(size=9),
        plot.title.position = "plot") +
  ggtitle("Basket, Basketless, All (n = 436 + 459)")+ 
  xlim(0, 15000) +
  ylim(0, 0.0012)

plt2 <- ggplot(res_basket_features %>%
                filter(gene %in% c(basket_nsb_up$gene, basket_nsb_down$gene)), 
              aes(x = width)) +
  geom_density(alpha = .5, aes(fill = change)) +
  geom_density(data = features, aes(x = width), linetype = "dashed", color ="grey") +
  theme_minimal() +
  theme(plot.title = element_text(size=9),
        plot.title.position = "plot") +
  ggtitle("Basket (n = 746 + 743)")+ 
  xlim(0, 15000) +
  ylim(0, 0.0012)

plt3 <- ggplot(res_basket_features %>%
                filter(gene %in% c(total_basket_nsb_up$gene, 
                                   total_basket_nsb_down$gene)), 
              aes(x = width)) +
  geom_density(alpha = .5, aes(fill = change)) +
  geom_density(data = features, aes(x = width), linetype = "dashed", color ="grey") +
  theme_minimal() +
  theme(plot.title = element_text(size=9),
        plot.title.position = "plot") +
  ggtitle("Basket, All (n = 194 + 172)")+ 
  xlim(0, 15000) +
  ylim(0, 0.0012)

ggarrange(plt2, plt3, plt1, legend = "bottom", common.legend = T, nrow = 1)
