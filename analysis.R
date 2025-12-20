#!/usr/bin/env Rscript

# Differential expression for GSE171546 RAW count files.
# 简单流程：找文件 -> 解析分组 -> 合并矩阵 -> DESeq2 对比 -> 导出 Excel 与火山图。

# --- 配置区域 ---------------------------------------------------------------
# 读取长行的 UTF-16LE 文件时，vroom 的默认缓冲区可能不足，先调大。
if (Sys.getenv("VROOM_CONNECTION_SIZE") == "") {
  Sys.setenv(VROOM_CONNECTION_SIZE = 5e6)  # 约 5 MB，避免 connection buffer 报错
}

# 必需 R 包
needed <- c("dplyr", "readr", "stringr", "purrr", "tibble", "DESeq2", "openxlsx", "ggplot2")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  stop("请先安装依赖: ", paste(missing, collapse = ", "))
}

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

# --- 1. 找到所有 Count 文件 --------------------------------------------------
count_dir <- "GSE171546_RAW"
files <- list.files(count_dir, pattern = "Count\\.txt\\.gz$", full.names = TRUE)
if (length(files) == 0) {
  stop("在 ", count_dir, " 未找到 *Count.txt.gz 文件")
}

# 从文件名提取样本名与分组： H_<sample>_Count.txt.gz
sample_info <- tibble(file = files) %>%
  mutate(
    sample = str_extract(basename(file), "(?<=H_).+(?=_Count\\.txt)"),
    replicate = str_extract(sample, "(?<=_)\\d+$"),
    condition = str_remove(sample, "_\\d+$")
  )

if (!any(sample_info$condition == "CON")) {
  stop("需要包含 CON 对照组")
}

# --- 2. 读取并合并计数矩阵 --------------------------------------------------
read_one <- function(path) {
  encodings <- c("UTF-16LE", "UTF-8", "latin1")
  ok <- NULL
  for (enc in encodings) {
    candidate <- read_tsv(path, locale = locale(encoding = enc), col_types = cols())
    if (ncol(candidate) >= 3) {
      ok <- candidate
      break
    }
  }
  if (is.null(ok)) {
    stop("文件无法正确解析: ", basename(path), "，请检查编码/分隔符")
  }

  sample_name <- str_extract(basename(path), "(?<=H_).+(?=_Count\\.txt)")
  tibble(gene_id = ok$gene_id, gene_name = ok$gene_name, !!sample_name := ok[[ncol(ok)]])
}

message("读取计数文件...")
count_list <- purrr::map(files, read_one)
counts <- purrr::reduce(count_list, dplyr::full_join, by = c("gene_id", "gene_name"))

# 重新排列表头，方便和样本信息一致
counts <- counts %>% select(gene_id, gene_name, all_of(sample_info$sample))

# 记录基因注释，避免后续对象被覆盖
gene_annotation <- counts %>% select(gene_id, gene_name) %>% distinct()

# --- 3. 构建 DESeq2 数据对象 -------------------------------------------------
coldata <- sample_info %>%
  transmute(
    sample = sample,
    condition = factor(condition, levels = c("CON", setdiff(unique(condition), "CON")))
  )
rownames(coldata) <- coldata$sample

count_matrix <- counts %>% select(-gene_name)
rownames(count_matrix) <- count_matrix$gene_id
count_matrix <- count_matrix %>% select(-gene_id)

# 过滤低表达：至少两个样本中计数 >= 10
keep <- rowSums(count_matrix >= 10) >= 2
message("基因数：", nrow(count_matrix), " -> 过滤后：", sum(keep))
count_matrix <- count_matrix[keep, ]
gene_annotation <- gene_annotation %>% filter(gene_id %in% rownames(count_matrix))

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(count_matrix)),
  colData = coldata,
  design = ~ condition
)

# 保留基因名用于输出
rowData(dds)$gene_name <- gene_annotation$gene_name[match(rownames(dds), gene_annotation$gene_id)]

message("运行 DESeq2...")
dds <- DESeq(dds)

# --- 4. 对比：每个实验组 vs CON --------------------------------------------
conditions <- setdiff(levels(coldata$condition), "CON")
if (length(conditions) == 0) stop("没有实验组可对比")

results_dir <- "results"
plot_dir <- file.path(results_dir, "volcano_plots")
dir.create(results_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)

wb <- createWorkbook()

for (grp in conditions) {
  message("对比 ", grp, " vs CON ...")
  res <- results(dds, contrast = c("condition", grp, "CON"))
  res_tbl <- as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(gene_name = rowData(dds)$gene_name[match(gene_id, rownames(dds))]) %>%
    relocate(gene_name, .after = gene_id) %>%
    arrange(padj, pvalue)

  sheet <- paste0(grp, "_vs_CON")
  addWorksheet(wb, sheet)
  writeData(wb, sheet, res_tbl)

  # 火山图（padj 和 |log2FC| 阈值为 0.05 / 1）
  plot_df <- res_tbl %>%
    mutate(
      padj_plot = pmax(padj, .Machine$double.xmin, na.rm = TRUE),
      neglog = -log10(padj_plot),
      sig = case_when(
        padj < 0.05 & abs(log2FoldChange) >= 1 ~ "padj<0.05 & |log2FC|>=1",
        padj < 0.05 ~ "padj<0.05",
        abs(log2FoldChange) >= 1 ~ "|log2FC|>=1",
        TRUE ~ "not sig"
      )
    )

  p <- ggplot(plot_df, aes(x = log2FoldChange, y = neglog, color = sig)) +
    geom_point(alpha = 0.7, size = 1.2) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "grey50") +
    scale_color_manual(values = c(
      "padj<0.05 & |log2FC|>=1" = "firebrick",
      "padj<0.05" = "steelblue",
      "|log2FC|>=1" = "darkorange",
      "not sig" = "grey70"
    )) +
    labs(title = paste0(grp, " vs CON"), x = "log2FC", y = "-log10(padj)", color = "标记") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(plot_dir, paste0(sheet, "_volcano.png")), p, width = 7, height = 6, dpi = 300)
}

# --- 5. 保存结果 ------------------------------------------------------------
excel_path <- file.path(results_dir, "deseq2_contrasts.xlsx")
saveWorkbook(wb, excel_path, overwrite = TRUE)

message("完成！结果文件：")
message("- ", excel_path)
message("- 火山图：", plot_dir)
