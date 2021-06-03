# this code generates the Nr5a1 TRAP-seq enrichment results

library(tidyverse)
library(ggrepel)

files <- list.files(path = "counts", full.names = TRUE)
samples <- str_extract(files, "Sample_[0-9]+")
metadata <- read_csv("metadata.csv") %>%
  mutate(Cells = factor(Cells, levels = c("Sup", "Bead"))) %>%
  mutate(Pair = str_extract(Sample_Name, "[0-9]$"))

# - Counts --------------------------------------------------------------------
counts <- map_dfc(files, read_tsv, skip = 4, col_names = FALSE) %>%
  select(`X1...1`, starts_with("X2")) %>%
  set_names("gene_id", samples)

genes <- read_csv("genes.csv") %>%
  filter(!duplicated(gene_id)) %>%
  select(gene_id, gene_name)


# - Positive controls ---------------------------------------------------------
counts %>%
  mutate_if(is.numeric, ~ .x / sum(.x) * 10^6) %>%
  filter(gene_id %in% c("ENSMUSG00000026751", "Cre", "Gfp_L10a")) %>%
  pivot_longer(-gene_id, names_to = "Sample_ID", values_to = "Counts") %>%
  mutate(gene_id = factor(
    gene_id,
    levels = c("ENSMUSG00000026751", "Cre", "Gfp_L10a"),
    labels = c("Nr5a1", "Cre", "Gfp_L10a")
  )) %>%
  left_join(select(metadata, Sample_ID, Cells, Pair), by = "Sample_ID") %>%
  ggplot(aes(x = Cells, y = Counts, group = Pair)) +
  geom_point(show.legend = FALSE, color = "black") +
  geom_line(show.legend = FALSE, color = "black") +
  scale_y_continuous(expand = c(0,0), trans = "log2") +
  theme_bw() +
  labs(y = "Counts per million", x = element_blank()) +
  facet_wrap(~gene_id) +
  theme(panel.grid = element_blank())


# - Heatmap -------------------------------------------------------------------
expression <- counts %>%
  mutate_if(is.numeric, ~ .x / sum(.x) * 10^6) %>%
  mutate_if(is.numeric, ~ .x > 1) %>%
  pivot_longer(-gene_id, names_to = "sample", values_to = "expressed") %>%
  group_by(gene_id) %>%
  summarize(in3 = sum(expressed) >= 3)

annotation <- data.frame(Cells = metadata$Cells)
rownames(annotation) <- metadata$Sample_ID

annotation_colors <- list(Cells = c(Bead = "#7fc97f", Sup = "gray"))

counts %>%
  mutate_if(is.numeric, ~ log2(.x / sum(.x) * 10^6 + 1)) %>%
  filter(gene_id %in% filter(expression, in3)$gene_id) %>%
  as.data.frame() %>%
  column_to_rownames("gene_id") %>%
  pheatmap::pheatmap(
    show_rownames = FALSE,
    annotation_col = annotation,
    annotation_colors = annotation_colors
  )


# - DESeq2 --------------------------------------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts %>% as.data.frame() %>% column_to_rownames("gene_id"),
  colData = metadata,
  design = ~ Pair + Cells
) %>%
  DESeq2::DESeq()

enrichment <- dds %>%
  DESeq2::results() %>%
  as.data.frame() %>%
  rownames_to_column("gene_id") %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

fpm <- dds %>%
  DESeq2::fpm() %>%
  as.data.frame() %>%
  rownames_to_column("gene_id")

fpm_means <- map(
  levels(metadata$Cells),
  ~ rowMeans(select(fpm, filter(metadata, Cells == .x)$Sample_ID))
) %>%
  as.data.frame() %>%
  set_names(levels(metadata$Cells)) %>%
  mutate(gene_id = fpm$gene)

# MA plot of enrichment by bead expression
enrichment %>%
  select(gene_id, log2FoldChange, padj) %>%
  left_join(fpm_means, by = "gene_id") %>%
  ggplot(aes(x = Bead, y = log2FoldChange, color = padj < 0.05)) +
  geom_point(show.legend = FALSE, alpha = 0.4) +
  scale_color_manual(values = c("black", "red")) +
  scale_x_continuous(trans = "log2", expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Bead expression (FPM)", y = "Fold enrichment (log2 bead/sup)") +
  theme_bw() +
  theme(panel.grid = element_blank())


# - Enrichment results --------------------------------------------------------
left_join(fpm_means, enrichment, by = "gene_id") %>%
  full_join(genes, ., by = "gene_id") %>%
  select(gene_name, gene_id, Bead, Sup, log2FoldChange, padj) %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
  rename("Enrichment" = log2FoldChange, "P" = padj) %>%
  write.csv(., "results.csv", row.names = FALSE, na = "")
