library(tidyverse)
library(ggrepel)

# read in count data
read_counts <- function(count_files) {
  samples <- str_extract(count_files, "Sample_[0-9]+")
  df <- map(count_files, ~read_tsv(.x, skip = 4, col_names = FALSE))
  df <- bind_cols(df) %>%
    select(`X1...1`, starts_with("X2")) %>%
    set_names("gene_id", samples)
  
  # add gene names
  tx2gene <- read_csv("genes.csv")
  df <- left_join(df, select(tx2gene, gene_id, gene_name), by = "gene_id")
  
  # use Ensembl name if no gene_name
  df <- mutate(df, gene_name = ifelse(is.na(gene_name), gene_id, gene_name))
  df <- dplyr::select(df, -gene_id)
  df <- dplyr::rename(df, "gene" = gene_name)
  return(df)
}

# make counts into cpm
make_cpm <- function(counts, log2 = FALSE) {
  if (log2 == TRUE) {
    df <- mutate_if(counts, is.numeric, ~ log2(.x / sum(.x) * 10^6 + 1))
  } else {
    df <- mutate_if(counts, is.numeric, ~ .x / sum(.x) * 10^6)
  }
  return(df)
}

find_expressed <- function(counts, threshold = 1, samples = 3) {
  cpm <- make_cpm(counts, log2 = FALSE)
  above_threshold <- mutate_if(cpm, is.numeric, ~ .x >= threshold)
  total_samples <- rowSums(select(above_threshold, -gene))
  return(cpm$gene[total_samples >= samples])
}

# plot library size for each sample
plot_library_size <- function(counts) {
  size <- colSums(select(counts, -gene))
  size <- as.data.frame(size) %>%
    rownames_to_column("sample") %>%
    mutate(size = size / 10^6)
  
  plt <- ggplot(size, aes(x = sample, y = size)) +
    geom_col() +
    geom_text(aes(label = round(size, 0), y = size + 1)) +
    labs(y = "Library size (M)", x = element_blank()) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  return(plt)
}

enrichment <- function(counts, metadata, genes, remove = NULL) {
  if (!is.null(remove)) {
    metadata <- filter(metadata, !(Sample %in% remove))
  }
  metadata <- select(metadata, Sample, Cells)
  metadata <- mutate(metadata, pair = rep(seq(nrow(metadata)/2), each = 2))
  
  # convert counts to cpm
  counts <- make_cpm(counts, log2 = TRUE)
  counts <- filter(counts, gene %in% genes)
  counts <- gather(counts, -gene, key = "Sample", value = "cpm")
  
  # add metadata to counts
  counts <- left_join(counts, metadata, by = "Sample")
  
  # generate enrichment score for each sample
  counts %>%
    group_by(pair, gene) %>%
    mutate("Enrichment" = cpm - lead(cpm)) %>%
    filter(Cells == "Bead") %>%
    ungroup() %>%
    select(-Cells, -pair, -cpm) %>%
    spread(key = "gene", value = "Enrichment") %>%
    as.data.frame()
}

# enrichment plot
plot_enrichment <- function(counts, metadata, genes, 
                            log2 = TRUE, pair = FALSE,
                            remove = NULL) {
  if (!is.null(remove)) {
    metadata <- filter(metadata, !(Sample %in% remove))
  }
  
  metadata <- select(metadata, Sample, Cells, Treatment)

  if (pair == TRUE) {
    metadata <- mutate(metadata, pair = rep(seq(nrow(metadata)/2), each = 2))
  }
  
  cpm <- make_cpm(counts, log2 = FALSE)
  cpm <- filter(cpm, gene %in% genes)
  cpm <- gather(cpm, -gene, key = "Sample", value = "cpm")
  cpm <- inner_join(cpm, metadata, by = "Sample")
  cpm <- mutate(cpm, Cells = factor(Cells, levels = c("Sup", "Bead")))
  cpm <- mutate(cpm, gene = factor(gene, levels = genes))
  
  
  # plot
  plt <- ggplot(cpm, aes(x = Cells, y = cpm, color = Treatment)) +
    facet_wrap(~gene) +
    theme_classic() +
    labs(y = "Counts per million (log2)", x = element_blank())
  
  if (log2 == TRUE) {
    plt <- plt + scale_y_continuous(trans = "log2")
  }
  if (pair == TRUE) {
    plt <- plt + geom_point() + geom_line(aes(group = pair))
  } else {
    plt <- plt + geom_jitter(height = 0)
  }
  return(plt)
}


# heatmap plot with dendrogram
plot_heatmap <- function(counts, metadata, genes, annotation = NULL,
                         show_genes = FALSE) {
  cpm <- make_cpm(counts, log2 = TRUE)
  cpm <- dplyr::filter(cpm, gene %in% genes & !duplicated(gene) & !is.na(gene)) %>% 
    as.data.frame()
  cpm <- column_to_rownames(cpm, "gene")

  # make annotation data.frame
  if (is.null(annotation)) {
    annotation_df <- data.frame(
      "Run" = metadata$Run,
      "Cells" = metadata$Cells,
      "Treatment" = metadata$Treatment,
      "Genotype" = metadata$Mutant
    )
  } else {
    annotation_df <- metadata[ , annotation] %>% as.data.frame()
  }
  rownames(annotation_df) <- metadata$Sample

  # plot
  pheatmap::pheatmap(cpm,
                     annotation_col = annotation_df,
                     show_rownames = show_genes)
}

# correlation plot with dendrogram
plot_correlation <- function(counts, metadata, genes, annotation = NULL) {
  cpm <- make_cpm(counts)
  cpm <- mutate_if(cpm, is.numeric, ~ log2(.x + 1))
  cpm <- filter(cpm, gene %in% find_expressed(counts) & !duplicated(gene) & !is.na(gene)) %>% 
    as.data.frame()
  cpm <- column_to_rownames(cpm, "gene")
  
  # calculate pairwise correlation
  search_grid <- expand.grid(colnames(cpm), colnames(cpm))
  results <- map(1:nrow(search_grid), ~ cor(cpm[, search_grid$Var1[.x]], 
                                            cpm[, search_grid$Var2[.x]]))
  search_grid$cor <- unlist(results)
  
  # reshape for plotting
  search_grid <- spread(search_grid, key = "Var2", value = "cor") %>%
    as.data.frame() %>%
    column_to_rownames("Var1")
  
  # make annotation data.frame
  if (is.null(annotation)) {
    annotation_df <- data.frame(
      "Run" = metadata$Run,
      "Cells" = metadata$Cells,
      "Treatment" = metadata$Treatment,
      "Genotype" = metadata$Mutant
    )
  } else {
    annotation_df <- metadata[ , annotation] %>% as.data.frame()
  }
  rownames(annotation_df) <- metadata$Sample
  
  # plot
  pheatmap::pheatmap(search_grid,
                     annotation_col = annotation_df,
                     annotation_row = annotation_df,
                     color = viridis::viridis(20))
}

# remove duplicate rownames by keeping highest expressing version of gene name
dedup_genes <- function(counts) {
  expr <- data.frame("rownumber" = seq(nrow(counts)), 
                     "total" = rowSums(select(counts, -gene)),
                     "gene" = counts$gene)
  expr <- arrange(expr, desc(total)) %>% filter(!duplicated(gene))
  expr <- arrange(expr, rownumber)
  return(counts[expr$rownumber, ])
}

# run DESeq2 with specified contrasts
DESeq_results <- function(dds, input, option) {
  if (option == "name") {
    df <- 
      dds %>%
      results(name = input) %>%
      as.data.frame() %>%
      rownames_to_column("gene")
  } else if (option == "contrast") {
    df <- 
      dds %>%
      results(contrast = input) %>%
      as.data.frame() %>%
      rownames_to_column("gene")
  }
  return(df)
}


# volcano plot
plot_volcano <- function(deseq_results, deseq_contrasts = NULL, genes = NULL,
                         n_col = NULL, label_genes = NULL) {
  if (!is.null(deseq_contrasts)) {
    deseq_results <- filter(deseq_results, contrast %in% deseq_contrasts) %>%
      mutate(contrast = factor(contrast, levels = deseq_contrasts))
  }
  # only keep selected genes
  if (!is.null(genes)) {
    deseq_results <- filter(deseq_results, gene %in% genes)
  } 
  
  # set axis limits
  min_pval <- min(deseq_results$padj, na.rm = TRUE)
  y_lim <- c(min_pval, 0)
  
  max_fc <- deseq_results %>%
    filter(padj < 0.05) %>% 
    summarize(max = max(abs(log2FoldChange), na.rm = TRUE)) %>%
    .$max
  x_lim <- c(-max_fc, max_fc)
  
  # plot
  plt <- ggplot(deseq_results, aes(x = log2FoldChange, y = -log10(padj),
                                   color = padj < 0.05)) +
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 0)) +
    geom_hline(aes(yintercept = -log10(0.05)), linetype = "dashed") +
    geom_point(show.legend = FALSE) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_continuous(limits = x_lim) +
    scale_color_manual(values = c("gray50", "red")) +
    theme_classic() +
    facet_wrap(~contrast, ncol = n_col) +
    labs(x = expression("Fold change (log"[2]*")"),
         y = expression("-log"[10]*italic(" P")*" value"))
  
  if (!is.null(label_genes)) {
    to_label <- filter(deseq_results, gene %in% label_genes)
    plt <- plt + 
      geom_text_repel(data = to_label, aes(label = gene), color = "black")
  }
  return(plt)
}

# correlation plot of de genes
plot_similarity <- function(deseq_results, genes, gene_names = FALSE,
                            show_legend = FALSE) {
  df <- filter(deseq_results, gene %in% genes)
  df <- select(df, gene, log2FoldChange, contrast)
  
  # make 0 the midpoint by chopping values at extremes
  center_fc <- function(values) {
    max_value <- max(values)
    min_value <- min(values)
    if (max_value > abs(min_value)) {
      new_values <- ifelse(values > abs(min_value), abs(min_value), values)
    } else {
      new_values <- ifelse(values < -max_value, -max_value, values)
    }
    return(new_values)
  }
  
  df$log2FoldChange <- center_fc(df$log2FoldChange)
  
  # set up data for heatmap
  df <- spread(df, key = "contrast", value = "log2FoldChange")
  df <- as.data.frame(df) %>% column_to_rownames("gene")
  
  # make red --> white --> blue color palette
  endcolors <- c("firebrick3", "white", "dodgerblue3")
  cols <- c(colorRampPalette(c(endcolors[1], endcolors[2]))(50), 
            colorRampPalette(c(endcolors[2], endcolors[3]))(51)[-1])
  
  pheatmap::pheatmap(df, show_rownames = gene_names, color = cols,
                     legend = show_legend,
                     border_color = NA)
}

# plot pairs
plot_pair <- function(deseq_results, contrast1, contrast2, genes = NULL,
                      xlab = NULL, ylab = NULL, center = FALSE,
                      limits = NULL, smooth = TRUE) {
  
  if (is.null(genes)) {
    genes <- filter(deseq_results, !is.na(padj)) %>%
      filter(!duplicated(gene)) %>% .$gene
  }
  # grab data
  df <- dplyr::filter(deseq_results, gene %in% genes)
  df <- select(df, gene, log2FoldChange, contrast)
  df <- spread(df, key = "contrast", value = "log2FoldChange")
  
  
  contrast1 <- enquo(contrast1)
  contrast2 <- enquo(contrast2)
  
  # plot
  plt <- ggplot(df, aes(x = !!contrast1, y = !!contrast2)) +
    geom_hline(aes(yintercept = 0)) +
    geom_vline(aes(xintercept = 0)) +
    geom_hline(aes(yintercept = log2(1.5)), linetype = "dashed") +
    geom_hline(aes(yintercept = -log2(1.5)), linetype = "dashed") +
    geom_vline(aes(xintercept = log2(1.5)), linetype = "dashed") +
    geom_vline(aes(xintercept = -log2(1.5)), linetype = "dashed") +
    geom_point() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_bw() +
    theme(panel.grid = element_blank())
  if (!is.null(xlab)) {
    plt <- plt + xlab(xlab)
  }
  if (!is.null(ylab)) {
    plt <- plt + ylab(ylab)
  }
  if (!is.null(limits)) {
    plt <- plt + 
      scale_x_continuous(expand = c(0, 0), limits = limits) + 
      scale_y_continuous(expand = c(0, 0), limits = limits)
  } else if (center == TRUE) {
    xlim <- c(-max(abs(select(df, !!contrast1))), max(abs(select(df, !!contrast1))))
    ylim <- c(-max(abs(select(df, !!contrast2))), max(abs(select(df, !!contrast2))))
    plt <- plt + 
      scale_x_continuous(expand = c(0, 0), limits = xlim) + 
      scale_y_continuous(expand = c(0, 0), limits = ylim)
  }
  if (smooth == TRUE) {
    plt <- plt + geom_smooth(method = "lm", fullrange = TRUE, se = FALSE, 
                color = "dodgerblue3")
  }
  return(plt)
}

plot_fchistogram <- function(comp, genes) {
  df <- filter(deseq_results, contrast == comp & gene %in% genes)
  
  limits <- c(-max(abs(df$log2FoldChange), na.rm = TRUE),
              max(abs(df$log2FoldChange), na.rm = TRUE))
  plt <- ggplot(df, aes(x = log2FoldChange)) +
    geom_histogram(fill = "black") +
    xlim(limits) +
    scale_y_continuous(expand = c(0, 0)) +
    xlab(NULL) + ylab(NULL) +
    theme_bw() +
    theme(panel.grid = element_blank())
  return(plt)
}

# correlation plot with dendrogram
plot_deseq_corr <- function(deseq_results, genes) {
  df <- filter(deseq_results, gene %in% genes)
  df <- select(df, gene, log2FoldChange, contrast)
  df <- spread(df, key = "contrast", value = "log2FoldChange")
  df <- as.data.frame(df) %>% column_to_rownames("gene")
  
  # calculate pairwise correlation
  search_grid <- expand.grid(colnames(df), colnames(df))
  results <- map(1:nrow(search_grid), ~ cor(df[, search_grid$Var1[.x]], 
                                            df[, search_grid$Var2[.x]]))
  search_grid$cor <- unlist(results)
  
  # reshape for plotting
  search_grid <- spread(search_grid, key = "Var2", value = "cor") %>%
    as.data.frame() %>%
    column_to_rownames("Var1")
  
  pheatmap::pheatmap(search_grid,
                     color = viridis::viridis(20))
}

# correlation test
test_correlation <- function(deseq_results, contrast1, contrast2, 
                             genes = NULL) {
  
  if (is.null(genes)) {
    genes <- filter(deseq_results, !is.na(padj)) %>%
      filter(!duplicated(gene)) %>% .$gene
  }
  # grab data
  df <- dplyr::filter(deseq_results, gene %in% genes)
  df <- select(df, gene, log2FoldChange, contrast)
  df <- spread(df, key = "contrast", value = "log2FoldChange")
  
  tst <- cor.test(df[, contrast1], df[, contrast2], method = "pearson")
  return(paste("r:", round(tst$estimate, 3), "| p:", round(tst$p.value, 4)))
}

# plot treatment
plot_treatment <- function(genes, tx = NULL, geno = NULL, order_tx = FALSE, 
                           order_genes = FALSE, ncol = NULL, colors = NULL,
                           mutant = FALSE) {
  fpm <- filter(fpm, gene %in% genes)
  fpm <- left_join(fpm, runs, by = "Sample")
  fpm <- filter(fpm, Cells == "Bead")
  
  if (!is.null(tx)) {
    fpm <- filter(fpm, Treatment %in% tx)
  }
  if (!is.null(geno)) {
    fpm <- filter(fpm, Mutant %in% geno)
  }
  if (order_tx == TRUE) {
    fpm <- mutate(fpm, Treatment = factor(Treatment, levels = tx))
  }
  if (order_genes == TRUE) {
    fpm <- mutate(fpm, gene = factor(gene, levels = genes))
  }
  
  # plot
  plt <- ggplot(fpm, aes(x = Treatment, y = fpm, fill = Treatment)) +
    geom_boxplot(show.legend = FALSE) +
    scale_y_continuous(trans = "log2") +
    facet_wrap(~gene, scales = "free_y", ncol = ncol) +
    labs(y = "Expression (FPM)", x = element_blank()) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 11))
  if (mutant == TRUE) {
    plt <- plt + facet_wrap(~ Mutant + gene, scales = "free_y", ncol = ncol)
  } 
  if (!is.null(colors)) {
    plt <- plt + scale_fill_manual(values = colors)
  }
  return(plt)
}


fc_heatmap <- function(deseq_results, contrasts, genes, gene_size = 12,
                       fill_limit = NULL) {
  if (is.null(fill_limit)) {
    fill_limit <- max(abs(filter(deseq_results, padj < 0.05)$log2FoldChange))
  }
  
  df <- deseq_results %>%
    filter(contrast %in% contrasts) %>%
    filter(gene %in% genes) %>%
    select(gene, log2FoldChange, contrast)
  df <- spread(df, key = "contrast", value = "log2FoldChange") %>%
    as.data.frame() %>%
    column_to_rownames("gene")
  
  # find distance
  gene_dist <- dist(df)
  gene_clusters <- hclust(gene_dist)
  gene_order <- gene_clusters$labels[gene_clusters$order]
  
  contrast_dist <- dist(t(df))
  contrast_clusters <- hclust(contrast_dist)
  contrast_order <- contrast_clusters$labels[contrast_clusters$order]
  
  df <- df %>% rownames_to_column("gene") %>%
    gather(-gene, key = "contrast", value = "fc") %>%
    mutate(contrast = factor(contrast, levels = contrast_order),
           gene = factor(gene, levels = gene_order))
  
  # plot
  plt <- ggplot(df, aes(x = contrast, y = gene, fill = fc)) +
    geom_tile() +
    scale_fill_gradient2(low = "firebrick3", mid = "white", high = "dodgerblue3",
                         name = "Fold Change\n(log2)",
                         limits = c(-fill_limit, fill_limit)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme_bw() +
    xlab(NULL) + ylab(NULL) +
    theme(axis.text = element_text(color = "black", size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size = gene_size),
          panel.grid = element_blank())
  return(plt)
}


gene_direction <- function(genes, comps = NULL, only_sig = FALSE) {
  if (!is.null(comps)) {
    deseq_results <- filter(deseq_results, contrast %in% comps)
  } else {
    comps <- names(contrasts)
  }
  
  deseq_results <- filter(deseq_results, gene %in% genes)
  output_names <- paste(rep(unique(deseq_results$contrast), each = 2), 
                        c("up", "down"))
  
  reg_genes <- function(name) {
    group <- str_extract(name, "^[\\S]+")
    if (str_detect(name, "up")) {
      filter(deseq_results, contrast == group & log2FoldChange > 0)$gene
    } else if (str_detect(name, "down")) {
      filter(deseq_results, contrast == group & log2FoldChange < 0)$gene
    }
  }
  result <- map(output_names, reg_genes)
  names(result) <- output_names
  return(result)
}

# pairplot of all correlations
plot_snspairplot <- function(genes, contrasts = NULL) {
  if (is.null(contrasts)) {
    to_plot <- unique(deseq_results$contrast)
    num <- length(to_plot)
  }
  
  df <- expand.grid(to_plot, to_plot) %>% as.data.frame()
  
  plots <- list()
  count <- 1
  for (i in 1:nrow(df)) {
    if (df[i, ]$Var1 != df[i, ]$Var2) {
      plots[[count]] <- ggplot(filter(deseq_results, contrast == df$Var1),
                               aes(x = log2FoldChange)) + geom_histogram() + theme_classic()
    } else {
      plots[[count]] <- 
        ggplot(filter(deseq_results, contrast == df$Var1), 
               aes(x = log2FoldChange)) + 
        geom_histogram() + theme_classic()
    }
    count <- count + 1
  }
  return(plots)
  #cowplot::plot_grid(plots)
}


venn_data <- function(up_down_list, direction) {
  keep <- up_down_list[which(str_detect(names(up_down_list), direction))]
  genes <- unique(unlist(keep))
  groups <- str_extract(names(keep), "[\\S]+")
  names(keep) <- groups
  
  remove_dups <- function(combo_df) {
    combo_df <- apply(combo_df, 1, sort) %>% t() %>% as.data.frame()
    combo_df <- combo_df[!duplicated(combo_df), ]
    return(combo_df)
  }
  
  intersection <- function(comp, n) {
    if (n == 2) {
      hits <- intersect(keep[[comp$V1, comp$V2]])
    } else if (n == 3) {
      hits <- intersect(keep[[comp$V1]], 
                        intersect(keep[[comp$V2]], keep[[comp$V3]])
      )
    } else if (n == 4) {
      hits <- intersect(keep[[comp$V1]], 
                       intersect(keep[[comp$V2]], 
                                 intersect(keep[[comp$V3]], keep[[comp$V4]])
                                )
      )
    }
  }
  
  # results
  if (length(groups) == 2) {
    combos <- expand.grid(groups, groups)
    combos <- remove_dups(combos)
  } else if (length(groups) == 3) {
    combos <- expand.grid(groups, groups, groups)
    combos <- remove_dups(combos)
  } else if (length(groups) == 4) {
    combos <- expand.grid(groups, groups, groups, groups)
    combos <- remove_dups(combos)
  }
  
  # sort combos to test broadest first
  uniques <- apply(combos, 1, function(x) length(unique(x))) %>% as.data.frame
  uniques <- rownames_to_column(uniques, "rownum")
  uniques <- arrange(uniques, desc(`.`))
  combos <- combos[uniques$rownum, ]
  
  # calculate the length of the gene intersection list
  intersection_length <- vector()
  for (i in 1:nrow(combos)) {
    hit_genes <- intersection(combos[i, ], ncol(combos))
    intersection_length[i] <- length(hit_genes[hit_genes %in% genes])
    genes <- genes[!(genes %in% hit_genes)]
  }
  
  # make a named vector for output
  mask <- apply(combos, 1, function(x) !duplicated(x)) %>% t()
  list_names <- vector()
  for (i in 1:nrow(combos)) {
    output <- vector()
    for (j in 1:ncol(combos)) {
      if (mask[i, j] == TRUE) {
        output <- c(output, as.character(combos[i, j]))
        }
      }
    list_names[i] <- paste(output, collapse = "&")
    }
  
  names(intersection_length) <- list_names
  intersection_length <- intersection_length[!duplicated(names(intersection_length))]
  return(intersection_length)
}


# - Gene ontology ------------------------------------------------------------
genes_to_csv <- function(comparison, folder) {
  genes <- filter(deseq_results, contrast == comparison) %>%
    filter(padj < 0.05) %>%
    .$gene
  if (length(genes) > 0) {
    print(paste('Writing CSV of', comparison))
    write.csv(genes, file = paste0(folder, comparison, "_de_genes.csv"),
              row.names = FALSE)
  }
}

read_go_results <- function(folder, df = FALSE) {
  files <- list.files(path = folder, pattern = "*results.txt", full.names = TRUE)
  samples <- str_extract(files, "([A-Za-z0-9]+)(?=_results)")
  column_names <- c("GO", "total", "hits", "expected",
                    "over/under", "enrichment", "P", "FDR")
  results <- map(files, ~ read_tsv(.x, skip = 11, col_types = c("cdddcddd")))
  results <- map(results, ~ set_names(., column_names))
  names(results) <- samples
  
  if (df == TRUE) {
    results <- map(1:length(results), ~ mutate(results[[.x]], 
                                               "contrast" = samples[.x])) %>%
      bind_rows()
  }
  return(results)
}

make_go_list <- function(go_terms) {
  print('Converting to GO list')
  terms <- unique(go_terms$GO)
  result <- map(terms, ~ filter(go_terms, GO == .x)$gene)
  result <- map(result, ~ unique(.x))
  names(result) <- terms
  return(result)
}

go_test <- function(comparison) {
  print(paste('Testing', comparison))
  deseq_results <- deseq_results %>%
    filter(padj < 0.05) %>%
    filter(contrast == comparison)
  genes <- deseq_results$padj
  names(genes) <- deseq_results$gene
  print(paste(length(genes), "DE genes"))
  
  if (length(genes) > 0) {
    # set up go data
    go_data <- new("topGOdata",
                   ontology = "BP",
                   allGenes = genes,
                   GO2genes = go_list,
                   annot = annFUN.GO2genes,
                   geneSel <- function(x) {x < 0.05},
                   nodeSize = 5)
    
    # run test
    results_ks <- runTest(go_data, algorithm = "classic", statistic = "ks")
    full_table <- GenTable(go_data, KS = results_ks, orderBy = "KS", 
                           topNodes = length(results_ks@score))
    full_table <- mutate(full_table, "contrast" = comparison)
    
    return(full_table)
  } else {
    return(NULL)
  }
}

geometric_mean <- function(nums) {
  n = length(nums)
  return(prod(nums)^(1/n))
}

harmonic_mean <- function(nums) {
  return(1/mean(1/nums))
}
