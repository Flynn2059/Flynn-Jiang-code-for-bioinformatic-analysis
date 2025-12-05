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
# 常规处理，barcode前加样本名避免重复，也便于后续分析设定分组信息
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
dds=DESeqDataSetFromMatrix(
  countData = pseudobulk_sample,
  colData   = metadata,
  design    = ~ group
)
# 过滤极低表达量的基因
dds=dds[rowSums(counts(dds)) > 10, ]
# 差异表达分析
dds=DESeq(dds)
# 提取 LPS vs CRLPS 的结果
res=results(dds, contrast = c("group", "LPS", "CRLPS"))
# 排序、查看显著基因
res=res[order(res$padj), ]
res=as.data.frame(res)
res_filtered=res %>%
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
vsd=vst(dds, blind = FALSE)
top_genes=rownames(res)[1:20]
mat=assay(vsd)[top_genes, ]
mat_z=t(scale(t(mat)))
# 设置热图颜色
col_fun <- colorRamp2(c(-1,0,2), c("#0f86a9","white","#FC8452"), transparency =0.2) # 有很多可选的热图配色，这里是其中的一个

# 自定义圆角图形函数
round_cell <- function(j,i,x,y,w,h,fill){
  grid.roundrect(x,y,w,h,r = unit(0.3,'snpc'),
                 gp = gpar(col = "white",fill = fill,lwd = 2))
}
# 设置热图颜色
col_fun <- colorRamp2(c(-1, 0, 2), c("#A5CC26", "white", "#FF7BAC"), transparency = 0.2)
col_fun <- colorRamp2(c(-1, 0, 2), c("#B3A9EB", "white", "#ffa500"), transparency = 0.2)

# 绘制热图
ht2 <- Heatmap(
   mat_z,
  name = "Expression",
  col = col_fun,
  
  # 设置格子样式
  cell_fun = round_cell,
  rect_gp = gpar(type = "none"),
  # 设置格子尺寸
  width = ncol( mat_z) * unit(0.6, "cm"),
  height = nrow( mat_z) * unit(0.6, "cm"),
  
  # 热图外观设置
  row_names_side="left",
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  
  # 使用圆角矩形（通过 raster 参数实现）
  raster_quality = 2,
  
  # 聚类设置（可选）
  cluster_rows = T,
  cluster_columns = F,
  
  # 颜色条设置
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 7)
  )
)
# 绘制热图
draw(ht2, padding = unit(c(2, 2, 2, 2), "mm"))

# 创建基因分组（每3个一组，共9组）
gene_groups <- c(rep(1:9, each = 3))
# group_colors <- sample(colors(), 9) 
# 自定义颜色
group_colors <- c("#D2749B","#F89C54","#F8AD20","#E8C559","#FDDB50","#E3E592","#93CC82","#9FD9DC","#A1D4F0","#C7C2E1")

# 创建顶部注释（基因分组）
top_annotation <- HeatmapAnnotation(
  cluster = anno_block(
    gp = gpar(fill = "white",col="white"),
    height = unit(0.25, "cm"),
    labels = row.names(data),
    labels_gp = gpar(col = "black", fontsize = 7)
  ),
  Group = anno_block(
    gp = gpar(fill = group_colors,col="white"),
    height = unit(0.1, "cm"),
    labels = NULL,
    labels_gp = gpar(col = "black", fontsize = 6)
  ),
  show_annotation_name = FALSE
)

# 设置热图颜色
col_fun <- colorRamp2(c(-1, 0, 2), c("#A5CC26", "white", "#FF7BAC"), transparency = 0.2)
col_fun <- colorRamp2(c(-1, 0, 2), c("#B3A9EB", "white", "#ffa500"), transparency = 0.2)

# 绘制热图
ht3 <- Heatmap(
   mat_z,
  name = "Expression",
  col = col_fun,
  
  # 设置格子样式
  cell_fun = round_cell,
  rect_gp = gpar(type = "none"),
  # 设置格子尺寸
  width = ncol( mat_z) * unit(0.6, "cm"),
  height = nrow( mat_z) * unit(0.6, "cm"),
  
  # 热图外观设置
  row_names_side="left",
  row_names_gp = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 45,
  
  # 使用圆角矩形（通过 raster 参数实现）
  raster_quality = 2,
  
  # 聚类设置（可选）
  cluster_rows = F,
  cluster_columns = F,
  
  # 列分割（用于分组显示）
  column_split = gene_groups,
  column_title = NULL, 
  
  # 顶部注释
  top_annotation = top_annotation,
  
  # 颜色条设置
  heatmap_legend_param = list(
    title = "Z-score",
    title_gp = gpar(fontsize = 9),
    labels_gp = gpar(fontsize = 7)
  )
)
# 绘制热图
draw(ht3, padding = unit(c(2, 2, 2, 2), "mm"))
