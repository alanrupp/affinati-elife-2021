library(tidyverse)
library(ggrepel)
library(wesanderson)

files <- list.files(pattern = "ReadsPerGene.out.tab", recursive = TRUE,
                    full.names = TRUE)
samples <- str_extract(files, "Sample_[0-9]{6}")

info <- read_csv("metadata.csv") %>%
  filter(!duplicated(Sample_ID)) %>%
  mutate(Cells = factor(Cells, levels = c("Sup", "Bead"))) %>%
  mutate(pair = rep(c(1,2), each = 2)) %>%
  mutate(pair = factor(pair))

# - Counts --------------------------------------------------------------------
counts <- map_dfc(files, read_tsv, skip = 4, col_names = FALSE) %>%
  select(`X1...1`, starts_with("X2")) %>%
  set_names("gene", samples)

genes <- read_csv("genes.csv") %>%
  filter(!duplicated(gene_name)) %>%
  select(starts_with("gene")) %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name))
genes[genes$gene_name == "Flpo", "gene_biotype"] <- "protein_coding"

counts <- counts %>%
  left_join(select(genes, gene_id, gene_name), by = c("gene" = "gene_id")) %>%
  mutate(gene_name = ifelse(is.na(gene_name), gene, gene_name)) %>%
  select(-gene) %>%
  rename(gene = gene_name)

# - Positive controls ---------------------------------------------------------
counts %>%
  mutate_if(is.numeric, ~ .x / sum(.x) * 10^6) %>%
  filter(gene %in% c("Lepr", "Cre", "Slc17a6", "Flpo", "Gfp_L10a")) %>%
  pivot_longer(-gene, names_to = "Sample_ID", values_to = "counts") %>%
  mutate(gene = factor(gene, levels = c("Lepr", "Cre", "Slc17a6",
                                  "Flpo", "Gfp_L10a"))) %>%
  left_join(select(info, Sample_ID, Run_ID, Cells, pair), by = "Sample_ID") %>%
  mutate(Cells = factor(Cells, levels = c("Sup", "Bead"))) %>%
  ggplot(aes(x = Cells, y = counts, group = pair, color = Run_ID)) +
  geom_point(show.legend = FALSE) +
  geom_line(show.legend = FALSE) +
  scale_y_continuous(trans = "log2") +
  scale_color_manual(values = wes_palette("Darjeeling1")) +
  theme_bw() +
  labs(y = "Counts per million", x = element_blank()) +
  facet_grid(~gene) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# - Heatmap ----------------------------------
library(pheatmap)
expression <-
  counts %>%
  mutate_if(is.numeric, ~ .x / sum(.x) * 10^6) %>%
  mutate_if(is.numeric, ~ .x > 1) %>%
  gather(-gene, key = "sample", value = "expressed") %>%
  group_by(gene) %>%
  summarize(in3 = sum(expressed) >= 3)

annotation <- data.frame(
  Run = info$Run_ID,
  Cells = info$Cells
)
rownames(annotation) <- info$Sample_ID

annotation_colors <- list(
  Cells = c(Bead = "#7fc97f", Sup = "gray"),
  Run = c("Run_2347" = wes_palette("GrandBudapest2")[1],
          "Run_2519" = wes_palette("GrandBudapest2")[2])
)

counts %>%
  mutate_if(is.numeric, ~ .x / sum(.x) * 10^6) %>%
  mutate_if(is.numeric, ~ log2(.x + 1)) %>%
  filter(gene %in% filter(expression, in3)$gene) %>%
  as.data.frame() %>%
  column_to_rownames("gene") %>%
  pheatmap(.,
           show_rownames = FALSE,
           annotation_col = annotation,
           annotation_colors = annotation_colors,
           cutree_cols = 2)

# - DESEq2 -------------------------------------------------------------------
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = counts %>% as.data.frame() %>% column_to_rownames("gene"),
  colData = info,
  design = ~ pair + Cells
) %>%
  DESeq2::DESeq()

enrich <- dds %>%
  DESeq2::results() %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  arrange(desc(padj < 0.05), desc(log2FoldChange))

fpm <- dds %>%
  DESeq2::fpm() %>%
  as.data.frame() %>%
  rownames_to_column("gene")

fpm_means <-
  map(unique(info$Cells),
      ~rowMeans(select(fpm, filter(info, Cells == .x)$Sample_ID))) %>%
  as.data.frame() %>%
  set_names(unique(info$Cells)) %>%
  mutate(gene = fpm$gene)

# Scatter plot of enrichment by bead expression
enrich %>%
  filter(gene %in% filter(genes, gene_biotype == "protein_coding")$gene_name) %>%
  select(gene, log2FoldChange, padj) %>%
  left_join(., fpm_means, by = "gene") %>%
  ggplot(aes(x = Bead, y = log2FoldChange, color = padj < 0.05)) +
  geom_point(show.legend = FALSE, alpha = 0.4, stroke = 0) +
  scale_color_manual(values = c("black", "red")) +
  scale_x_continuous(trans = "log2", expand = c(0,0), limits = c(0.02, 4000)) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 4.3)) +
  geom_text_repel(data = data.frame(
    Bead = c(6.6, 0.7, 243.6, 56.5, 464.4),
    log2FoldChange = c(0.345, 0.352, 1.24, 1.20, 3.023),
    gene = c("Lepr", "Cre", "Slc17a6", "Flpo", "Gfp_L10a")),
    aes(label = gene), color = "black") +
  labs(x = "Bead expression (FPM)",
       y = "Fold enrichment (log2 bead/sup)",
       title = "Lepr-enriched genes") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# - Results ------------------------------------------------------------------
# format results for output
results <- left_join(enrich, fpm_means, by = "gene") %>%
  select(gene, Bead, Sup, log2FoldChange, padj) %>%
  rename(Gene = gene, Enrichment = log2FoldChange, P = padj)

write.csv(results, "results.csv", row.names = FALSE)
write.csv(fpm, "cpm.csv", row.names = FALSE)
