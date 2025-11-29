# ==== 一些我们喜欢的设置 ====
Sys.setenv(LANGUAGE="en")
options(StringAsFactors=FALSE)
rm(list=ls());gc()

setwd("~/BioInformatics Codes/")
getwd()

library(qs)
library(dplyr)
library(DESeq2)
library(Seurat)
library(ggpubr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(DoubletFinder)

# ==== 加载数据 ====
load("t_cell.RData")
print(t_cell)
table(t_cell@meta.data$celltype)
colnames(t_cell)=paste0(t_cell@meta.data$orig.ident,"_",colnames(t_cell))

# ==== FindMarkers示例 ====
# 这段代码会报错，看看就好
# DEGs=FindMarkers(t_cell,group.by = "Group",ident.1 = "LPS",ident.2 = "CRLPS",
#                 slot="counts",test.use="DESeq2")
counts=GetAssayData(t_cell,layer = "counts")
counts=as.data.frame(counts)
sample_ids=sub("_(.*)$", "", colnames(counts))
# 确保分组向量顺序与表达矩阵列顺序一致
names(sample_ids)=colnames(counts)
# 将 data.frame 转成数值矩阵
counts_mat=as.matrix(counts)
# 按样本对细胞进行求和：先转置再用 rowsum 按“细胞 → 样本”分组
pseudobulk_sample=t(rowsum(t(counts_mat), group = sample_ids))

# ==== 准备PseudoBulk ====
metadata=data.frame(
  row.names =unique(sample_ids),
  group=as.factor(c(rep("CRLPS",3),rep("LPS",3)))
)
dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData   = metadata,
  design    = ~ group
)
# 过滤极低表达量的基因
dds <- dds[rowSums(counts(dds)) > 10, ]
# 差异表达分析
dds <- DESeq(dds)

# 提取 LPS vs CRLPS 的结果
res <- results(dds, contrast = c("group", "LPS", "CRLPS"))

# 排序、查看显著基因
res <- res[order(res$padj), ]
res <- as.data.frame(res)

res_filtered <- res %>%
  filter(padj < 0.05 & abs(log2FoldChange) >= 0.8)

# --- EnhancedVolcano Plot ---
library(EnhancedVolcano)
library(Seurat)
plot=EnhancedVolcano(res,
                     lab = rownames(res),
                     x = 'log2FoldChange',
                     y = 'padj',
                     pCutoff = 0.05,
                     FCcutoff = 0.8,
                     title = '',subtitle = "")+NoGrid()
plot

# --- Heatmap of top 20 DEGs ---
library(pheatmap)
vsd <- vst(dds, blind = FALSE)
top_genes <- rownames(res)[1:20]
mat <- assay(vsd)[top_genes, ]
mat_z <- t(scale(t(mat)))
pheatmap(mat_z, 
         annotation_col = as.data.frame(colData(dds)[, "group", drop = FALSE]),
         scale = "column",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_names_row = FALSE,
         show_colnames = FALSE, 
         main = "Top 20 DEGs")


