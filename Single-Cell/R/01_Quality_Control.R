# ==== 加载配置和工作环境 ====
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
gc()
setwd("/Volumes/FlynnDisk/Red/Respond/")
getwd()
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
seurat_obj=ReadMtx(
  mtx="matrix.mtx.gz",
  cells = "barcodes.tsv.gz",
  features = "features.tsv.gz",
  feature.column = 1
)
seurat_obj=CreateSeuratObject(seurat_obj,project = "Respond",min.cells = 10,min.features = 50) 
print(seurat_obj)

# ==== 质量控制_去除低质量的细胞 ====
## 计算线粒体基因占比
seurat_obj[["percent_mt"]]=PercentageFeatureSet(seurat_obj,pattern = "^MT")
## 计算核糖体基因占比
seurat_obj[["percent_ribo"]]=PercentageFeatureSet(seurat_obj, pattern = "^RPS|^RPL")
# 计算红细胞基因占比
seurat_obj[["percent_RBC"]] <- PercentageFeatureSet(seurat_obj, pattern = "^HB[AB]")
# 绘制小提琴图
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt", "percent_RBC"), 
        pt.size = 0,ncol = 2)
# 进行质控，标准如下
seurat_obj=subset(seurat_obj,subset = nFeature_RNA>500 & percent_mt<15  & percent_RBC<1)
print(seurat_obj)
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent_mt"), 
        pt.size = 0,ncol = 3)

# 简单跑一下流程
seurat_obj=NormalizeData(seurat_obj,normalization.method = "LogNormalize",
                         scale.factor = median(seurat_obj@meta.data$nCount_RNA)) %>% FindVariableFeatures() %>% 
  ScaleData() %>% RunPCA() 
ElbowPlot(seurat_obj,ndims=50)
seurat_obj=RunUMAP(seurat_obj,dims=1:8)
# 找合适的pK数
possible_pK=paramSweep(seurat_obj,PCs=1:8,sct = F)
## PCs：主成分，上面的值
## sct：有没有使用SCTransform算法
summarize_pK=summarizeSweep(possible_pK,GT = F)
## GT：ground truth，有没有人工标记出单细胞，大多数数据都没有，所以F
best_pK=find.pK(summarize_pK)
best_pK
best_pK=0.2
# 标记双细胞并且查看结果
nExp=round(0.08*ncol(seurat_obj)) # 华大C4平台的估计双细胞数目
seurat_obj=DoubletFinder::doubletFinder(seurat_obj,PCs = 1:8,pN = 0.25,pK = best_pK,
                                        nExp = nExp,reuse.pANN = F,sct = F)
colnames(seurat_obj@meta.data)
table(seurat_obj@meta.data[["DF.classifications_0.25_0.2_1926"]])
seurat_obj=subset(seurat_obj,DF.classifications_0.25_0.2_1926=="Singlet")
print(seurat_obj)
qsave(seurat_obj,file = "Respond_post_QC.qs")

