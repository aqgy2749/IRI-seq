<style>
.jp-RenderedHTMLCommon, .jp-RenderedHTMLCommon p, .jp-RenderedHTMLCommon li,
.text_cell_render, .text_cell_render p, .text_cell_render li {
  font-size: 12px !important;
  line-height: 1.35 !important;
}
.jp-RenderedHTMLCommon h1, .text_cell_render h1 { font-size: 20px !important; }
.jp-RenderedHTMLCommon h2, .text_cell_render h2 { font-size: 16px !important; }
.jp-CodeCell pre, .input_area pre { font-size: 11px !important; }
@media print {
  body { font-size: 10pt !important; }
  h1 { font-size: 16pt !important; }
  h2 { font-size: 13pt !important; }
  pre, code { font-size: 8.5pt !important; }
}
</style>
# 将实验记录转为 HTML 并打印
1. 在终端进入本文件所在目录。
2. 运行：`jupyter nbconvert --to html EXPERIMENT_LOG.md`，或在 Jupyter 中打开后导出 HTML。
3. 下载 HTML 到本地并使用浏览器打开。
4. 按 `Ctrl+P`（macOS 使用 `Command+P`），选择打印机或“另存为 PDF”。

# IRI-seq 复现实验记录
- 项目：`IRI-seq`
- 记录范围：2026-07-03 至 2026-07-06
- 当前状态：FY1_3 raw FASTQ 到 cDNA 表达矩阵和 bead-bead connection 的复现流程已搭建；metadata 注释对齐存在样本字段与 barcode 前缀错位问题。
- 依据：本目录 notebook、辅助脚本、作者 GitHub 脚本、GEO processed 文件及本轮运行输出。
- 主要目录：
  - 代码目录：`/p1/zulab_users/jtian/my_jupyter/IRI-seq`
  - 作者脚本：`/p1/zulab_users/jtian/my_jupyter/IRI-seq/script_from_github`
  - FY1_3 输入：`/p2/zulab/jtian/data/IRISeq/demo_FY1_3/data`
  - FY1_3 输出：`/p2/zulab/jtian/data/IRISeq/demo_FY1_3/output/demo-FY1_3-pipeline`

## 2026-07-03

**内容：** Fig. 2C / GEO processed 数据初步复现

**提交：** 当前工作区修改，尚未提交

**目的：** 根据论文、GEO 页面和作者 GitHub 代码，判断 Fig. 2C 对应数据，并尝试用作者 processed expression matrix 复现正文中的两张图。

**修改：**
- 新建或更新 `IRIS/demo_f2_processed.ipynb`。
- 读取 GEO processed expression matrix：`GSE270383_count.mtx.gz`、`GSE270383_barcodes.tsv.gz`、`GSE270383_genes.tsv.gz`。
- 读取 `GSE270383_meta_data.csv.gz` 用于 annotation 和坐标信息。
- 按作者代码思路尝试重画 expression UMAP 和 spatial reconstruction 图。

**结果：**
- 确认 `GSE270383_count.mtx.gz` 是全数据集 processed matrix，不是单个样本矩阵。
- 发现只用 processed expression matrix 很难直接复现正文 Fig. 2C，因为正文图依赖作者后续整合、annotation 和空间重构坐标。
- 初步暴露出 `metadata` 中 barcode index 前缀和 `orig.ident` 样本标签可能不完全一致的问题。

## 2026-07-04

**内容：** FY1_3 raw FASTQ 到重构图的全流程 notebook

**提交：** 当前工作区修改，尚未提交

**目的：** 最大限度沿用作者 GitHub 逻辑，构建 FY1_3 / Section 2 从 FASTQ 到 cDNA 表达矩阵、bead-bead connection 和重构图的完整复现 notebook。

**修改：**
- 新建 `IRIS/demo-FY1_3-pipeline.ipynb`。
- 固定输入目录为 `demo_FY1_3/data`，输出目录为 `demo_FY1_3/output/demo-FY1_3-pipeline`。
- 建立 FASTQ 映射：
  - `SRR29481311_1/2.fastq.gz` -> `20231221_FY1_3_cDNA.R1/R2.fastq.gz`
  - `SRR29481287_1/2.fastq.gz` -> `20231221_FY1_3_connection.R1/R2.fastq.gz`
  - 将 `SRR29481287_2.fastq.gz` 同时作为作者脚本需要的 `R3` 读入。
- cDNA 分支按作者 EasySpatial 顺序执行：
  - spatial/UMI barcode attach
  - cutadapt 去 polyA
  - STAR 比对
  - samtools 过滤
  - barcode+UMI+位置去重
  - gene counting
  - `Generate_adata` 生成 `adata_full.h5ad`
- bead interaction 分支按作者脚本执行：
  - `UMI_barcode_extraction.py`
  - `spatial_barcode_extraction.py`
  - `Remove_duplicate_barcode.py`
- 增加 processed GEO 文件验证章节，用于和 raw FASTQ 复现结果对比。

**主要参数：**

```python
UMI_LIMIT_FOR_ADATA = 50
EXPRESSION_MIN_COUNTS = 400
HVG_N_TOP_GENES = 10000
EXPRESSION_PCS = 25
EXPRESSION_NEIGHBOR_PCS = 20
EXPRESSION_UMAP_MIN_DIST = 0.1
MIN_CONNECTION_UMIS = 8
RECONSTRUCTION_UMAP_N_NEIGHBORS = 19
RECONSTRUCTION_UMAP_MIN_DIST = 0.23
RECONSTRUCTION_UMAP_SPREAD = 0.5
RECONSTRUCTION_UMAP_EPOCHS = 500
```

**结果：**
- cDNA 分支生成 `cDNA/Adata/adata_full.h5ad`。
- `adata_full.h5ad` 中保留 `7902` 个 barcode、`28590` 个 gene。
- 表达 UMAP 进一步按 `total UMI >= 400` 过滤后剩 `5573` 个点。
- bead connection 分支生成 deduplicated spatial connection 文件：
  - `beads/Spatial_barcode_rmdup/20231221_FY1_3_connection.spatial.csv.gz`
- 作者 GEO processed connection 文件和复现 connection 文件行数均为 `31942236`。

**注意事项：**
- 作者 bead interaction 脚本存在 FASTQ header 换行导致 5-line pseudo-FASTQ 的兼容问题，notebook 中增加了修复逻辑。
- `seurat_v3` HVG 依赖 `scikit-misc`，当前环境缺失时会 fallback 到 Scanpy `seurat` flavor。

## 2026-07-05

**内容：** FY1_3 pipeline 各 cell 逻辑解释与中间统计核对

**提交：** 当前工作区修改，尚未提交

**目的：** 解释 notebook 各 cell 输出，明确 cDNA barcode、UMI、STAR 唯一比对、gene assignment、dedup 和 `adata_full` barcode 数的来源。

**修改：**
- 逐步解释 `demo-FY1_3-pipeline.ipynb` 中 cDNA pipeline、metadata loading、expression UMAP、bead pipeline、spatial reconstruction 和 validation 的代码逻辑。
- 增加或明确 Jupyter display 设置，避免 DataFrame 中路径和长 barcode 被省略。
- 读取并解释 STAR log：
  - `cDNA/Sam_STAR/20231221_FY1_3_cDNALog.final.out`

**关键结果：**
- barcode extraction：
  - raw reads：`402742745`
  - barcode/UMI attach 后 reads：`168997196`
- cutadapt 后 STAR input reads：`141470693`
- STAR uniquely mapped reads：`25180667`，比例 `17.80%`
- barcode+UMI+位置去重后 reads：`12992450`
- gene assigned reads：`10688260`
- `Generate_adata` 检测到 cell barcode：`215310`
- `UMI_count > 50` 后保留：`7902`
- 这些有效 barcode 捕获约 `96.93%` 的 UMI。

**解释结论：**
- `215310` 来自 gene count summary 中所有出现过的 `Cell_name` 数量，不是最终用于分析的高质量 bead/cell 数。
- `7902` 是 `Cell_name` 按 gene-level UMI 汇总后 `UMI_count > 50` 的 barcode 数。
- `5573` 是在 `7902` 基础上再执行 `EXPRESSION_MIN_COUNTS = 400` 后的表达 UMAP 输入点数。

## 2026-07-06

**内容：** metadata barcode 匹配策略调整、重绘结果图和 validation 更新

**提交：** 当前工作区修改，尚未提交

**目的：** 检查 `GSE270383_meta_data.csv.gz` 中 FY1_3 barcode/index 前缀与 `orig.ident` 的关系，解释 annotation 匹配少、重构图结构不对的问题。

**修改：**
- 更新 `IRIS/demo-FY1_3-pipeline.ipynb` 的 metadata loading cell：
  - 原逻辑优先使用 `orig.ident == FY1_3`
  - 调整为优先使用 `metadata.index` 包含 `FY1_3`
  - 打印 index 命中数、`orig.ident` 命中数和两者交集数
- 新增辅助 runner：
  - `IRIS/_rerun_fy1_3_from_metadata.py`
  - 用于复用已有中间结果，重跑 metadata、expression UMAP annotation、主图和 validation。
- 新增 validation 更新脚本：
  - `IRIS/_update_fy1_3_validation_index_metadata.py`
  - 用于刷新 index 匹配后的 validation summary。
- 在 `demo-FY1_3-pipeline.ipynb` 最后追加 metadata debug cell：
  - 读取 `GSE270383_meta_data.csv.gz`
  - 筛选 barcode/index 名字包含 `FY1_3` 的行
  - 打印 `orig.ident` 分布和 `Annotation` 分布
  - 显示 `barcode_index`、`orig.ident`、`Annotation` 和 UMAP 坐标列。

**关键统计：**

```text
metadata index contains FY1_3: 3920
metadata orig.ident == FY1_3: 3920
both: 123
```

**processed expression matrix 对比：**

```text
reproduced adata_full barcode count: 7902
processed matrix total barcode count: 128405
processed barcodes index contains FY1_3: 3920
metadata orig.ident == FY1_3: 3920
compact overlap reproduced vs processed all: 3934
compact overlap reproduced vs processed index contains FY1_3: 3920
compact overlap reproduced vs metadata orig.ident == FY1_3: 125
```

**重绘结果：**
- 表达图注释上 `3920 / 5573` 个点。
- 空间重构图注释上 `3918 / 10700` 个点。
- 重绘输出：
  - `figures/FY1_3_expression_umap_by_annotation.png`
  - `figures/FY1_3_spatial_reconstruction_by_annotation.png`
  - `figures/FY1_3_expression_and_spatial_two_panel.png`

**validation 结果：**

```text
reproduced_connection_unique_molecules: 31942236
processed_connection_unique_molecules: 31942236
reproduced_edges_min_umi: 225728
processed_edges_min_umi: 225728
interaction_edge_jaccard_min_umi: 1.0
interaction_edge_UMI_pearson: 1.0
interaction_edge_UMI_spearman: 1.0
cDNA_common_barcodes_for_total_UMI: 3920
cDNA_total_UMI_pearson: 0.9905197416474466
cDNA_total_UMI_spearman: 0.9952071822781279
```

**结论：**
- `SRR29481287` bead interaction FASTQ 到 deduplicated connection 文件的复现结果与作者 GEO processed connection 文件高度一致，行数和 edge-level UMI 统计均一致。
- cDNA 侧的主要问题不是 barcode normalization 本身，而是作者 processed matrix 的 barcode/index 前缀和 metadata 中 `orig.ident` 样本身份字段存在错位。
- 按 index 前缀 `FY1_3` 匹配可以获得较高 barcode overlap，但 annotation 中混入大量 `section1` 和 `section3_4`，生物学上不适合作为 FY1_3 / Section 2 注释。
- 按 `orig.ident == FY1_3` 选择 annotation 生物学上更合理，但和 raw FASTQ 复现出的 barcode 只 overlap 约 `125` 个。
- 因此目前不能同时获得高 barcode overlap 和可信 FY1_3 annotation，这也是重绘 expression UMAP 和 reconstruction 图缺少正文结构的主要原因之一。

## 当前问题与后续方向

**cDNA barcode 匹配问题：**
- `GSE270383_count.mtx.gz`、`GSE270383_barcodes.tsv.gz`、`GSE270383_genes.tsv.gz` 是全数据集 processed matrix。
- 其中 barcode/index 名字包含 `FY1_3` 的列有 `3920` 个，和 raw FASTQ 复现的 `7902` 个 barcode 中的 `3920` 个可以匹配。
- 但 metadata 中真正 `orig.ident == FY1_3` 的 `3920` 行，与 raw FASTQ 复现 barcode 只匹配约 `125` 个。
- 后续分析应将 `metadata_by_index` 和 `metadata_by_origident` 分开：
  - `metadata_by_index` 用于技术性 barcode overlap 检查。
  - `metadata_by_origident` 用于生物学 FY1_3 / Section 2 annotation 检查。

**reconstruction 图结构问题：**
- 作者 CPU reconstruction notebook 使用 `log10(n_umi) >= 0.9`、`n_neighbors=19`、`min_dist=0.23`。
- 作者 GPU reconstruction notebook 使用 `log10(n_umi) >= 0.2`、`n_neighbors=70`、`min_dist=0.65`。
- 当前 notebook 主要按 CPU 参数实现，并为可运行性使用 PCA/SVD 后再 UMAP；这与正文最终图可能使用的流程不完全一致。
- 后续需要单独比较 CPU/GPU reconstruction 参数、是否直接对 connection matrix UMAP、以及是否需要作者保存的最终 spatial coordinates。

**建议下一步：**
- 先用 processed metadata 中 `orig.ident == FY1_3` 的 `UMAP1_spatial/UMAP2_spatial` 直接画作者 processed 坐标，作为生物学正确的 FY1_3 / Section 2 reference。
- 再将 raw/processed connection 生成的 reconstruction UMAP 与该 reference 做 Procrustes 或 KNN preservation 对比。
- 保留当前 `demo-FY1_3-pipeline.ipynb` 作为 raw FASTQ 全流程复现实验本，不把 index 匹配结果误认为最终生物学 annotation。
