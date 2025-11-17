# ==== 加载配置和工作环境 ====
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls()); gc()
setwd("~/Flow_Cytometry/")
getwd()

suppressPackageStartupMessages({
  library(qs)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(SCP)
  library(Seurat)
  library(cowplot)
  library(ggplot2)
  library(ggpubr)
  library(patchwork)
  library(stringr)
  library(tidyverse)
  library(RColorBrewer)
})

# ==== 我们假设你已经完成了数据读入和处理 ====
# 我们需要读入的表达矩阵是每行一个基因，每列一个样本的表达矩阵
# 使用经过荧光补偿后的表达矩阵，不要做数据放缩/Scaling
# 这里我们假设处理好的表达矩阵是data，请再检查一次
head(data)

# ==== 当作Seurat对象读取进去 ====
seurat_obj=CreateSeuratObject(counts = data)
print(seurat_obj)
# 跳过质控部分，应该在FlowJo里完成
# 需要去除碎片，多胞和死细胞；如果有需要也可以先筛选出自己感兴趣的细胞类型
# 可以根据选择补充metadata信息，比如orig.ident
# seurat_obj@meta.data$orig.ident= #后面输入一个字符串向量，可以考虑使用paste0()根据样本信息生成
# ==== 降维聚类 ====
seurat_obj=NormalizeData(seurat_obj, normalization.method = 'CLR', margin = 2) 
VariableFeatures(seurat_obj) = rownames(seurat_obj[["RNA"]]) # 使用所有的标记进行分析
seurat_obj=ScaleData() 
seurat_obj=RunPCA()
print(seurat_obj)
# 一般来说不需要矫正批次效应，因为是同一个人同一台机器同一批样本染的，同样或大致相同的荧光补偿也应用于所有样本
seurat_obj=RunUMAP(seurat_obj,dims = 1:5) # 对于n个标的数据，dims应当输入1:n-1
seurat_obj=FindNeighbors(seurat_obj,dims = 1:5) # 对于n个标的数据，dims应当输入1:n-1
seurat_obj=FindClusters(seurat_obj,resolution = 0.05)
table(seurat_obj@meta.data$seurat_clusters)
# ====  可视化 ====
plot1=DimPlot(seurat_obj,group.by="seurat_clusters",reduction = "umap",label=T,raster = F)
# 为了更好的可视化效果，我们设定raster=F，如果需要加快出图速度可以设置为T

# ==== 细胞注释 ====
# 细胞注释的marker需要根据自己的流式染的marker去填写，这里只做示例
# 也可以不写，用rownames(seurat_obj)展示所有的marker
known_markers=list("T/NK"=c('CD3E','CD4','CD8A'),
                   B=c("CD79A","MS4A1"),
                   Plasma=c("MZB1"),
                   Epithelium=c('EPCAM',"KRT8","KRT18"),
                   Tumor=c("MKI67","FOXJ1","SOX2","SOX9"),
                   "Mono/Macro"=c("CD14","CD68","CD163"),
                   Mast=c("KIT"),
                   Fibroblast=c("COL5A2",'PDGFRB'),
                   Pericyte=c("CSPG4","RGS5"),
                   Endothelium=c("PECAM1")) 
# 经典气泡图，看不同cluster的marker gene表达情况
plot2=DotPlot(object = seurat_obj,
              features = known_markers, # 或者rownames(seurat_obj)展示所有的marker
              scale=T,
              group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
plot1|plot2
# 邱同学风格的细胞注释
meta_supp = data.frame(seurat_cluster = 0:(length(unique(seurat_obj$seurat_clusters)) - 1), celltype = NA)
meta_supp[meta_supp$seurat_cluster %in% c(0), 'celltype'] = 'B'
meta_supp[meta_supp$seurat_cluster %in% c(4), 'celltype'] = 'Plasma'
meta_supp[meta_supp$seurat_cluster %in% c(1), 'celltype'] = 'T/NK'
meta_supp[meta_supp$seurat_cluster %in% c(5), 'celltype'] = 'Stromal'
meta_supp[meta_supp$seurat_cluster %in% c(3), 'celltype'] = 'Mono/Macro'
meta_supp[meta_supp$seurat_cluster %in% c(2), 'celltype'] = 'Tumor'

for (i in 1:nrow(meta_supp)) {
  seurat_obj@meta.data[which(seurat_obj$seurat_clusters == meta_supp$seurat_cluster[i]), 'celltype_major'] = meta_supp$celltype[i]
}
Idents(seurat_obj) <- 'celltype_major'

# 看看注释情况
plot3=CellDimPlot(
  seurat_obj,
  group.by = "celltype_major",
  theme_use = "theme_blank",
  xlab = "UMAP1",
  ylab = "UMAP2",
  label = TRUE,           
  label_insitu = TRUE,
  show_stat = F,      
  legend.position = "none" 
)
plot4=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "celltype_major")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
plot3|plot4
table(seurat_obj@meta.data$celltype_major)

# ==== 不要忘记保存注释后的对象 ====
qsave(seurat_obj, file = 'flow_cytometry_data_seurat_obj_post_annotation.qs')
