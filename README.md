# project_for_chen

#现在我有这个文件夹，我希望用R语言实现以下功能，发顶刊：
#1.选出文件名中有“Count”的文件，他们是一个个样本的count表达矩阵
#2.现在，文件名中“H_”到“.txt”之间的部分是样本号，其中CON是对照组，其他是实验组。比如，“24H”是实验24h组。
#3.我需要：计算每个时间段实验组和对照组的差异基因，画好火山图，并导出excel表。要有logfc和矫正p值。

## 分析脚本
根目录下的 `analysis.R` 复现了以上需求，采用 DESeq2 构建严谨的差异表达流程：

1. 自动发现 `GSE171546_RAW` 目录中包含 `Count` 的 UTF-16LE 计数文件，并从文件名解析样本编号、分组和重复编号。
2. 合并基因计数矩阵，过滤低表达基因（至少两个样本中计数 ≥ 10），并以 CON 作为参考组拟合 DESeq2 模型。
3. 针对每个时间点（CLP24H、CLP48H、CLP72H）与对照组的对比，输出包含 log2FC、p 值和校正后 p 值的差异表达表格，并绘制火山图。
4. 所有结果保存在 `results/` 目录下：`deseq2_contrasts.xlsx`（每个对比一张工作表）以及 `results/volcano_plots/` 中的 PNG 图片。

### 环境准备
- 安装 R（建议 ≥ 4.3）。
- 安装所需 R 包：`dplyr`, `readr`, `stringr`, `purrr`, `tibble`, `DESeq2`, `openxlsx`, `ggplot2`。示例命令：
  ```r
  install.packages(c("dplyr", "readr", "stringr", "purrr", "tibble", "openxlsx", "ggplot2"))
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("DESeq2")
  ```
- 如读取 UTF-16LE 计数文件时报 `connection buffer` 不足，可在运行前设置更大的缓冲区（脚本中已默认提高，如需手动设置可执行）：  
  ```r
  Sys.setenv(VROOM_CONNECTION_SIZE = 5e6)  # 约 5 MB
  ```

### 运行
在仓库根目录执行：
```bash
Rscript analysis.R
```
脚本会输出进度信息，并在 `results/` 中生成 Excel 表和火山图。
