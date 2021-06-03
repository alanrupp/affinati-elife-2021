library(DESeq2)
library(tidyverse)
library(ggrepel)
source("functions.R")

# - Read in data -------------------------------------------------------------
runs <- read_csv("runs.csv")
count_files <- list.files("counts", full.names = TRUE)

runs <- mutate(runs, Treatment = factor(Treatment, levels = c("PBS", "Leptin")))

# read in count data
counts <- read_counts(count_files)

# - Basic quality control -----------------------------------------------------
# plot library size for each sample
plot_library_size(counts)

# positive control genes
plot_enrichment(counts, runs, c("Lepr", "Cre", "Gfp_L10a"), pair = TRUE)
enrichment(counts, runs, c("Lepr", "Cre", "Gfp_L10a")) %>%
  gather(-Sample, key = "gene", value = "enrichment") %>%
  ggplot(aes(x = Sample, y = enrichment)) +
  geom_col() +
  facet_wrap(~gene)

# clustered heatmap of log2 cpm for top expressed genes
plot_heatmap(counts, runs, find_expressed(counts))

# correlation matrix of log2 cpm for top expressed genes
plot_correlation(counts, runs, find_expressed(counts))


# - Run DESeq2 to find DE genes -----------------------------------------------
enrichment_info <- runs %>%
  mutate(pair = as.character(rep(seq(n()/2), each = 2)))

# remove duplicate gene names
counts <- dedup_genes(counts)

dds <- DESeqDataSetFromMatrix(
  countData = counts %>% select(gene, enrichment_info$Sample) %>%
    as.data.frame() %>% column_to_rownames("gene"),
  colData = enrichment_info,
  design = ~ pair + Cells
) %>%
  DESeq()

enrichment <- results(dds, contrast = c("Cells", "Bead", "Sup")) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

ggplot(enrichment, aes(x = log2FoldChange, y = -log10(padj),
                       color = padj < 0.05)) +
  geom_hline(aes(yintercept = 0)) +
  geom_vline(aes(xintercept = 0)) +
  geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
  geom_point(show.legend = FALSE) +
  geom_label_repel(data = filter(enrichment, gene %in% c("Cre", "Gfp_L10a", "Lepr")),
                   aes(label = gene), color = "black") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_color_manual(values = c("gray50", "red4")) +
  theme_classic() +
  labs(x = expression("Enrichment (log"[2]*" bead / sup)"),
       y = expression("-log"[10]*italic(" P")*" value"))

# get mean expression
cpm <- fpm(dds)
cpm_means <- map(
  unique(enrichment_info$Cells), 
  ~ rowMeans(cpm[, filter(enrichment_info, Cells == .x)$Sample])
) %>%
  bind_cols() %>%
  set_names(unique(enrichment_info$Cells)) %>%
  mutate("gene" = rownames(cpm))

# - Write results -------------------------------------------------------------
select(enrichment, gene, log2FoldChange, padj) %>%
  inner_join(cpm_means, by = "gene") %>%
  select(gene, Bead, Sup, log2FoldChange, padj) %>%
  dplyr::rename("Enrichment" = log2FoldChange, "P" = padj) %>%
  write.csv("results.csv", row.names = FALSE)
