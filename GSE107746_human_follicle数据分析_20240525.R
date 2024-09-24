# 聚类降维 ========================================================================
rm(list = ls())
gc()
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 3 * 1024^3)
setwd("/home/dengxy/scRNA-seq/代码复现/GSE107746_Folliculogenesis/")

library(Seurat)
library(dplyr)
follicle <- read.table("GSE107746_Folliculogenesis_FPKM.log2.txt", header = TRUE, row.names = 1,sep = "\t")
follicle <- CreateSeuratObject(counts = follicle)  # 26364x151
follicle[["percent.mt"]] <- PercentageFeatureSet(follicle, pattern = "^MT")
VlnPlot(follicle, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

follicle <- NormalizeData(follicle) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)  %>%
  ScaleData(VariableFeatures(follicle))  %>%
  RunPCA()
ElbowPlot(follicle)  # 10 PC
DimPlot(follicle,reduction = "pca") + NoLegend()

follicle <- FindNeighbors(follicle)
follicle <- FindClusters(follicle,resolution = 3)  # 0.01 分了2个亚群  # 2~3-OK,5以下
# follicle <- RunTSNE(follicle,dims = 1:10)
# DimPlot(follicle,reduction = "tsne")
follicle <- RunUMAP(follicle,dims = 1:10)
DimPlot(follicle,reduction = "umap",label = T)
# follicle@meta.data$RNA_snn_res.2.5 <- NULL

### 乔杰院士2018年的人卵母细胞-颗粒细胞数据
# Transcriptome Landscape of Human Folliculogenesis Reveals Oocyte and Granulosa Cell Interactions.
mk_oo <- c('ZP2','DDX4','SYCP3','SOX30',
           'ZAR1','DAZL','YBX2','LHX8',
           'NOL4','CTCTF',
           'GPD1','NTF4','LCP2')
DotPlot(follicle,features = mk_oo)+RotatedAxis()

mk_GC <- c('CYP11A1','STAR','INHBA','AMH',
           'ZEB2','CD44',
           'TST','RBM24')   # 'BNIPL','CDCA3',
DotPlot(follicle,features = mk_GC)+RotatedAxis()

DotPlot(follicle,features = c('ZP1','ZP2','ZP3','ZP4','FSHR','AMHR2'))+RotatedAxis()
VlnPlot(follicle,features = c('ZP1','ZP2','ZP3','ZP4','FSHR','AMHR2'))+RotatedAxis()

# I just feel something is wrong!