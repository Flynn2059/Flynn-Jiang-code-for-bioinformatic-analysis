# ==== 加载配置和工作环境 ====
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
gc()
setwd("/Volumes/FlynnDisk/Red/")
getwd()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(harmony)
library(patchwork)
library(RColorBrewer)

# ==== 读取数据 ====
non_respond=qread("./Non_Respond/Non_Respond_post_QC.qs")
respond=qread("./Respond/Respond_post_QC.qs")
seurat_list=list(non_respond,respond)
seurat_obj=Reduce(merge,seurat_list)
print(seurat_obj)

# ==== 降维聚类+去除批次效应 ====
seurat_obj[["RNA"]]=JoinLayers(seurat_obj[["RNA"]])
seurat_obj=NormalizeData(seurat_obj,normalization.method = "LogNormalize",
                       scale.factor = median(seurat_obj@meta.data$nCount_RNA))
seurat_obj=FindVariableFeatures(seurat_obj,nfeatures=3000)
# 我们发现在亚群注释的过程中，如果按照上面的结果直接去做细胞注释
# 会有一些奇奇怪怪的细胞高表达非经典神经内分泌marker
# 我们觉得这样的结果不反映我们的生物学现象，而且这些细胞恰好高表达核糖体基因
# 所以我们希望去除核糖体基因和细胞周期基因对降维聚类的影响
# 计算细胞周期评分
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
seurat_obj[["percent_ribo"]]=PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
# 回归掉不感兴趣的变量
seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score","percent_ribo"))
# 使用HVG去跑PCA
seurat_obj=RunPCA(seurat_obj)
ElbowPlot(seurat_obj,ndims = 50)
seurat_obj=RunHarmony(seurat_obj,group.by.vars="orig.ident")
seurat_obj=RunUMAP(seurat_obj,dims=1:8,verbose = T,reduction = "harmony")
seurat_obj=FindNeighbors(seurat_obj,dims = 1:8,reduction = "harmony")
seurat_obj=FindClusters(seurat_obj,resolution = 0.1)
table(seurat_obj@meta.data$seurat_clusters)
plot1=DimPlot(seurat_obj,reduction = "umap",group.by = "orig.ident",label = T)
plot2=DimPlot(seurat_obj,reduction = "umap",group.by = "seurat_clusters",label = T)+NoLegend()
plot1|plot2
# ==== 细胞注释 ====
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
plot3=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()


plot3=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "seurat_clusters")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
plot2|plot3
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
plot4=DimPlot(seurat_obj,group.by = "celltype_major",label = T)&NoLegend()
plot5=DotPlot(object = seurat_obj,
              features = known_markers,
              scale=T,
              group.by = "celltype_major")+
  scale_color_gradientn(colors=brewer.pal(9,"Blues"))+
  theme_pubr()+
  theme(axis.text.x = element_text(angle=90)) & NoLegend()
plot4|plot5
table(seurat_obj@meta.data$celltype_major)

# ==== 不要忘记保存注释后的对象 ====
qsave(seurat_obj, file = 'seurat_obj_post_annotation.qs')

# ==== 这里顺带把不同种类的细胞也subset了 ====
T_NK=subset(seurat_obj,subset=celltype_major=="T/NK")
qsave(T_NK,file = "T_NK.qs")
Mono_Macro=subset(seurat_obj,subset = celltype_major=="Mono/Macro")
qsave(Mono_Macro,file = "Mono_Macro.qs")

