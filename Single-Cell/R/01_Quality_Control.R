# ==== 加载配置和工作环境 ====
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls());gc()
setwd("/Volumes/FlynnDisk/")
getwd()
list.files()
library(qs)
library(dplyr)
library(Seurat)
library(ggpubr)
library(cowplot)
library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(DoubletFinder)
# ==== 加载数据 ====
seurat_obj=Read10X_h5("./GSM8288503/filtered_feature_bc_matrix.h5")
seurat_obj=CreateSeuratObject(seurat_obj,project = "GSM8288503",min.cells = 10,min.features = 50) 
colnames(seurat_obj)=paste0(seurat_obj@meta.data$orig.ident,"_",colnames(seurat_obj))
print(seurat_obj)

# ==== 质量控制_去除低质量的细胞 ====
## 计算线粒体基因占比
seurat_obj[["percent_mt"]]=PercentageFeatureSet(seurat_obj,pattern = "^MT")
## 计算核糖体基因占比
seurat_obj[["percent_ribo"]]=PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
# 计算红细胞基因占比
seurat_obj[["percent_RBC"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[AB]")
# 绘制频数分布图
hist(seurat_obj@meta.data$nCount_RNA,xlim = c(0,5000),breaks = 2000)
abline(v=700,col="red")
hist(seurat_obj@meta.data$nFeature_RNA,xlim=c(0,1000),breaks=100) # 每个细胞的RNA数
abline(v=400,col="red")
hist(seurat_obj@meta.data$percent_mt,xlim=c(0,100),breaks=10)
abline(v=20,col="red")# 线粒体基因占比

# 进行质控，标准如下
seurat_obj=subset(seurat_obj,subset = nCount_RNA>700 &
                    nFeature_RNA>400 &
                    percent_mt<50 &
                    percent_RBC<0.1)
print(seurat_obj)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt","percent_ribo"), 
        pt.size = 0,ncol = 2)

# 简单跑一下流程
seurat_obj=NormalizeData(seurat_obj,normalization.method = "LogNormalize",
                         scale.factor = median(seurat_obj@meta.data$nCount_RNA)) %>% FindVariableFeatures() %>% 
  ScaleData() %>% RunPCA() 
ElbowPlot(seurat_obj,ndims=50)
seurat_obj=RunUMAP(seurat_obj,dims=1:20)
seurat_obj=FindNeighbors(seurat_obj,dims = 1:20) %>% FindClusters(resolution=0.1)
# 找合适的pK数
possible_pK=paramSweep(seurat_obj,PCs=1:20,sct = F)
## PCs：主成分，上面的值
## sct：有没有使用SCTransform算法
summarize_pK=summarizeSweep(possible_pK,GT = F)
## GT：ground truth，有没有人工标记出单细胞，大多数数据都没有，所以F
best_pK=find.pK(summarize_pK)
best_pK
best_pK=0.03
# 标记双细胞并且查看结果
doublet_rate=(ncol(seurat_obj)/1000)*0.008
nExp=round(doublet_rate*ncol(seurat_obj)) # 10X平台的估计双细胞数目
seurat_obj=DoubletFinder::doubletFinder(seurat_obj,PCs = 1:20,pN = 0.25,pK = best_pK,
                                        reuse.pANN = F,nExp = nExp,sct = F)
# 自动找到DoubletFinder的列
doublet_info = grep("^DF\\.classifications_", 
                colnames(seurat_obj@meta.data), 
                value = TRUE)
table(seurat_obj@meta.data[[doublet_info]])
# 拿出Singlet
singlet = rownames(seurat_obj@meta.data)[
  seurat_obj@meta.data[[doublet_info]] == "Singlet"
]
seurat_obj <- subset(seurat_obj, cells = singlet)
print(seurat_obj)
qsave(seurat_obj,file = "GSM8288503_post_QC.qs")
