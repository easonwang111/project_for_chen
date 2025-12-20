#!/usr/bin/env Rscript

# Differential expression analysis pipeline for GSE171546 RAW counts.
# Steps implemented per README requirements:
# 1. Discover Count files and parse sample metadata from the filename.
# 2. Build a combined count matrix.
# 3. Run DESeq2 contrasts comparing each time point to controls.
# 4. Export per-contrast tables (with log2FC and adjusted p values) to Excel.
# 5. Generate volcano plots for each contrast.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(DESeq2)
  library(openxlsx)
  library(ggplot2)
})

# Validate package availability early for reproducibility.
required_packages <- c("dplyr", "readr", "stringr", "purrr", "tibble", "DESeq2", "openxlsx", "ggplot2")
missing_pkgs <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing required packages: ", paste(missing_pkgs, collapse = ", "), ". Install them before running this script.")
}

# Paths and discovery
count_dir <- "GSE171546_RAW"
count_files <- list.files(count_dir, pattern = "Count\\.txt\\.gz$", full.names = TRUE)
if (length(count_files) == 0) {
  stop("No count files found in ", count_dir, ". Ensure the raw files are present and match the expected pattern.")
}

# Extract sample metadata from filenames.
sample_info <- tibble(
  file = count_files,
  sample = str_extract(basename(count_files), "(?<=H_).+(?=_Count\\.txt)")
) %>%
  mutate(
    replicate = str_extract(sample, "(?<=_)\\d+$"),
    condition = str_remove(sample, "_\\d+$"),
    condition = factor(condition, levels = unique(condition)),
    sample = factor(sample, levels = sample)
  )

if (!any(sample_info$condition == "CON")) {
  stop("Control samples (CON) are required but were not found in the discovered files.")
}

# Utility: read a single count file using the correct encoding (UTF-16LE).
read_count_file <- function(path) {
  df <- read_tsv(
    path,
    locale = locale(encoding = "UTF-16LE"),
    col_types = cols()
  )
  sample_name <- str_extract(basename(path), "(?<=H_).+(?=_Count\\.txt)")
  tibble(
    gene_id = df$gene_id,
    gene_name = df$gene_name,
    !!sample_name := df[[ncol(df)]]
  )
}

# Build the combined count table (wide format).
message("Reading count files and building expression matrix ...")
count_tables <- map(sample_info$file, read_count_file)
counts_merged <- reduce(count_tables, full_join, by = c("gene_id", "gene_name")) %>%
  distinct(gene_id, .keep_all = TRUE)

# Order columns to match sample metadata.
counts_merged <- counts_merged %>% select(gene_id, gene_name, all_of(as.character(sample_info$sample)))

# Prepare DESeq2 dataset.
gene_annotation <- counts_merged %>% select(gene_id, gene_name) %>% distinct()
count_matrix <- counts_merged %>% select(-gene_name)
rownames(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix %>% select(-gene_id)

coldata <- sample_info %>%
  select(sample, condition) %>%
  mutate(
    sample = as.character(sample),
    condition = factor(condition, levels = c("CON", setdiff(levels(condition), "CON"))),
    condition = relevel(condition, ref = "CON")
  )
rownames(coldata) <- coldata$sample

# Filter low-expressed genes: keep genes with at least 10 counts in at least two samples.
keep_genes <- rowSums(count_matrix >= 10) >= 2
message("Genes before filtering: ", nrow(count_matrix), "; after filtering: ", sum(keep_genes))
count_matrix <- count_matrix[keep_genes, ]
gene_annotation <- gene_annotation %>% filter(gene_id %in% rownames(count_matrix))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_matrix)),
  colData = coldata,
  design = ~ condition
)

# Retain gene names for reporting.
rowData(dds)$gene_name <- gene_annotation$gene_name[match(rownames(dds), gene_annotation$gene_id)]

message("Running DESeq2 ...")
dds <- DESeq(dds)

contrast_levels <- setdiff(levels(coldata$condition), "CON")
if (length(contrast_levels) == 0) {
  stop("No experimental conditions found beyond controls; nothing to contrast.")
}

results_dir <- file.path("results")
plot_dir <- file.path(results_dir, "volcano_plots")
dir.create(results_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

# Helper to format and annotate DESeq2 results.
extract_results <- function(dds_obj, group) {
  contrast <- c("condition", group, "CON")
  res <- lfcShrink(dds_obj, contrast = contrast, type = "normal")
  res_df <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(
      gene_name = rowData(dds_obj)$gene_name[match(gene_id, rownames(dds_obj))]
    ) %>%
    relocate(gene_name, .after = gene_id) %>%
    arrange(is.na(padj), padj, pvalue)
  res_df
}

# Volcano plot generator.
plot_volcano <- function(res_df, group) {
  plot_df <- res_df %>%
    mutate(
      plot_padj = ifelse(is.na(padj), NA_real_, pmax(padj, .Machine$double.xmin)),
      neg_log10_padj = -log10(plot_padj),
      significance = case_when(
        padj < 0.05 & abs(log2FoldChange) >= 1 ~ "padj < 0.05 & |log2FC| ≥ 1",
        padj < 0.05 ~ "padj < 0.05",
        abs(log2FoldChange) >= 1 ~ "|log2FC| ≥ 1",
        TRUE ~ "not significant"
      )
    )

  ggplot(plot_df, aes(x = log2FoldChange, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.7, size = 1) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") +
    scale_color_manual(values = c(
      "padj < 0.05 & |log2FC| ≥ 1" = "firebrick",
      "padj < 0.05" = "steelblue",
      "|log2FC| ≥ 1" = "darkorange",
      "not significant" = "grey70"
    )) +
    labs(
      title = paste0(group, " vs CON"),
      x = "log2 fold change",
      y = "-log10(adjusted p-value)",
      color = "Significance"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

message("Computing contrasts and exporting results ...")
wb <- createWorkbook()

for (grp in contrast_levels) {
  grp_label <- paste0(grp, "_vs_CON")
  res_tbl <- extract_results(dds, grp)

  # Add to workbook
  addWorksheet(wb, grp_label)
  writeData(wb, sheet = grp_label, x = res_tbl)

  # Volcano plot
  volcano <- plot_volcano(res_tbl, grp)
  ggsave(
    filename = file.path(plot_dir, paste0(grp_label, "_volcano.png")),
    plot = volcano,
    width = 7,
    height = 6,
    dpi = 300
  )
}

excel_out <- file.path(results_dir, "deseq2_contrasts.xlsx")
saveWorkbook(wb, excel_out, overwrite = TRUE)

message("Analysis complete. Results saved to:")
message("- ", excel_out)
message("- Volcano plots in ", plot_dir)
