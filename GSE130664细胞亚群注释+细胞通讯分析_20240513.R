# GSE130664
# Ovary 聚类注释 OK -----
rm(list=ls())
gc()
setwd("/home/dengxy/scRNA-seq/代码复现/GSE130664-Cynomolgus Monkey Ovary/")
# install.packages('Seurat')
library(Seurat)
library(dplyr)  
# read.delim（），用于读取以制表符（tab）分隔的文本文件的函数。该函数用于从文本文件中读取数据并创建一个数据框（data frame）
GSE130664_barcode_information <- read.delim("GSE130664_barcode_information.txt")
GSE130664_merge_UMI_count <- read.delim("GSE130664_merge_UMI_count.txt", row.names=1)
#View(GSE130664_barcode_information)
#colnames(GSE130664_merge_UMI_count)
ovary <- CreateSeuratObject(counts = GSE130664_merge_UMI_count, project = "ovary", min.cells = 3, min.features =200)
ovary[["percent.mito"]] <- PercentageFeatureSet(object = ovary, pattern = "^MT")

# ovary@assays$RNA@data[grep("MT",rownames(ovary@assays$RNA@data)),]%>%rownames()  # 从ovary对象的RNA数据中筛选出行名包含字符串"MT"的行，并返回这些行的行名

# ovary[["RNA"]]$counts[grep("MT",rownames(ovary@assays$RNA@data)),]%>%rownames()
# LayerData(ovary, assay="RNA", layer='counts')
# GetAssayData(ovary, assay="RNA", slot='data')

VlnPlot(ovary, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 
head(ovary[["percent.mito"]])
head(ovary@meta.data, 5)
#ovary.mt=ovary[["percent.mito"]]
plot1 <- FeatureScatter(ovary, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(ovary, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2  #CombinePlots(plots = list(plot1, plot2))

# Filter cells
# 创建的 Ovary 对象为 19008x3369
# 原文保留了 2,601 single cells (418 oocytes and 2,183 somatic cells)
# ovary[["RNA"]]@features %>% rownames()  # 查看基因名称  19008
# ovary[["RNA"]]@features %>% rownames() %>% head(10)
# ovary[["RNA"]]@cells %>% rownames()     # 查看细胞数量  3369
table(ovary$orig.ident)
# OF1 OF2 OF3 OF4 YF1 YF2 YF3 YF4 
# 449 416 423 575 361 325 370 450

mean(ovary$nFeature_RNA) # 2793.654
mean(ovary$percent.mito) # 0.2062707  # 计算nCount_RNA的平均数
median(ovary$nCount_RNA) # 39970      # 计算nCount_RNA中位数
table(ovary$nFeature_RNA > 700 )  # nFeature_RNA < 700有386个
table(ovary$nCount_RNA > 300000)    # nCount_RNA < 3000有130个
table(ovary$percent.mito > 1) # percent.mito > 1有39个

VlnPlot(ovary,features = 'nCount_RNA', y.max = 530000)
VlnPlot(ovary,features = 'nFeature_RNA')
VlnPlot(ovary,features = 'percent.mito', y.max = 1)

table(ovary$nFeature_RNA > 700 &
        ovary$nFeature_RNA < 8000 &
        ovary$nCount_RNA > 3000 &
        ovary$nCount_RNA < 300000 &
        ovary$percent.mito < 1) # 2602
table(ovary$nFeature_RNA > 700 &
        ovary$nCount_RNA > 3000 ) # 2958
table(ovary$nFeature_RNA > 700 &
        #ovary$nFeature_RNA < 8000 &
        ovary$nCount_RNA > 3000 &
        ovary$nCount_RNA < 530000 &
        ovary$percent.mito < 1)  # 2867

# ovary <- subset(ovary, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & nCount_RNA < 500000) #涂师姐的参数
# ovary <- subset(ovary, subset = nFeature_RNA > 200 & nCount_RNA > 200 & nCount_RNA < 500000 & ovary.mt < 1)  
# ovary <- subset(ovary, subset = nFeature_RNA > 700 & nCount_RNA > 3000) #参考原文献
# ovary <- subset(ovary, subset = nFeature_RNA > 700 & nFeature_RNA < 8000 & nCount_RNA > 3000 & nCount_RNA < 300000 & percent.mito < 1)   # 19008x2602
# ovary <- subset(ovary, subset = nFeature_RNA > 700 & nFeature_RNA < 12000 & nCount_RNA > 3000 &nCount_RNA < 530000 & percent.mito < 0.8)  # 19008x2867
ovary <- subset(ovary, subset = nFeature_RNA > 700 & nCount_RNA > 3000 & nCount_RNA < 530000 & percent.mito < 1)  # 19008x2829
VlnPlot(ovary, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 

# Normalize the data 
ovary <- NormalizeData(object = ovary, normalization.method = "LogNormalize", scale.factor = 10000)
# 归一化后的数据存在data里

# 原文保留了 2,601 single cells (418 oocytes)
mk_oo <- c('LMOD3','RBM46','NETO1',  # novel gene identified by this article
           "SYCP3","DDX4","GDF9","ZP3", # mark gene
           "NLRP14","NLRP5","MYOCOS","NPM2","ZAR1") # selected this time
# 提取表达至少一个基因的细胞，可以使用|（或）操作符 
ovary_oo <- subset(ovary, subset = LMOD3 >1 | RBM46 >1 | NETO1 >1 |
                     SYCP3 > 1 | DDX4 > 1 | GDF9 > 1 | ZP3 >1 |
                     NLRP14 >1 | NLRP5 >1 | MYOCOS > 1| NPM2 > 1| ZAR1 > 1)  # 19008x413;  不对ovary进行subset时，有19008x471
VlnPlot(ovary_oo, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 
# 结合ovary_oo反推质控标准
table(ovary_oo$nFeature_RNA > 700 & # 不变
        #ovary_oo$nFeature_RNA < 12000 &
        ovary_oo$nCount_RNA > 3000 & # 不变
        ovary_oo$nCount_RNA < 530000 & 
        ovary_oo$percent.mito < 1) # 413
table(ovary$nFeature_RNA > 700 & # 不变
        #ovary$nFeature_RNA < 12000 &
        ovary$nCount_RNA > 3000 & # 不变
        ovary$nCount_RNA < 530000 & 
        ovary$percent.mito < 1) # 2867

# Find top 2000 variable genes
ovary <- FindVariableFeatures(ovary, selection.method = "vst", nfeatures = 2000)
# Plot variable genes
top10 <- head(VariableFeatures(ovary), 10)  # top 10 variable genes
plot1 <- VariableFeaturePlot(ovary)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1+plot2  #CombinePlots(plots = list(plot1, plot2))

# Scale data based on all genes
ovary <- ScaleData(ovary, features = rownames(ovary))
# Use ID'd 2000 Variable genes to run PCA analysis
ovary <- RunPCA(ovary, features = VariableFeatures(object = ovary))
DoHeatmap(ovary,features = VariableFeatures(ovary) %>% head(20))
# Plot heatmaps based on PCs
DimHeatmap(ovary, dims = 1:10, cells = 500, balanced = TRUE)
# Elbow plot and jackstraw plots are used to determine the number of PC to use 
ElbowPlot(ovary)  


#设置了并行计算的策略为多进程，并指定了使用 4 个进程。。
#install.packages("future")
#install.packages("doParallel")
library(future)
library(doParallel)
plan("multisession", workers = 4)   # future::plan("multisession", workers = 4)
# ovary <- JackStraw(ovary, num.replicate = 100)
# ovary <- ScoreJackStraw(ovary, dims = 1:20)
# JackStrawPlot(ovary, dims = 1:20)

# 原文分了 7 个细胞类群
ovary <- FindNeighbors(ovary, dims = 1:7)
ovary <- FindClusters(ovary, resolution = 0.17)
# 原文筛选了 418 oocytes
table(ovary$seurat_clusters)  #查看各聚类的细胞数量
#prop.table(table(ovary$seurat_clusters)) #查看各聚类的细胞比例
round(prop.table(table(ovary$seurat_clusters)) * 100, 3)

# As for oocyte
rm(mk_oo)
mk_oo <- c('LMOD3','RBM46','NETO1',  # novel gene identified by this article
           "SYCP3","DDX4","GDF9","ZP3", # mark gene
           "NLRP14","NLRP5","MYOCOS","NPM2","ZAR1") # selected this time
DotPlot(ovary,features = mk_oo)+RotatedAxis()  # 全都OK，特异性很好
# 食蟹猴 ovary Marker
rm(mk)
mk <- c("TCF21","COL1A2","STAR",  # SC
        "AMH","WT1","INHA","CYP19A1", # GCs
        "CD3D","KLRB1",  # NKT
        "SYCP3","DDX4","GDF9","ZP3",# oocyte
        "CDH5","VWF",  # EC
        "DES","ACTA2",  # SMCs
        "CD68","CD14")  # Macrophages

DotPlot(ovary,features = mk)+RotatedAxis()
#               0   1   2   3   4   5   6  ：找 oocyte 接近418的
# 1:7+0.17    1146  487  449  366  209  118   92 OK:oo-366
# 1:8+0.15    1153  466  462  367  208  119   92 OK:oo-367
# 1:9+0.13    1048  506  452  367  208  194   92 OK:oo-367
# 1:10+0.08   1270  494  367  248  210  192   86 OK:oo-367
# 1:12+0.08   992   559  531  367  269   123  26 
# 1:14+0.08   996   561  444  367  269   204  26
# 1:15+0.07   1559  441  366  272  120   83   26
# 1:17+0.06   1574  437  366  269  112   83   26 x
# 1:20+0.055  1496  444  367  295  109   83   73 x

ovary <- RunTSNE(ovary, dims = 1:9)  #1:8
DimPlot(ovary, reduction = "tsne", pt.size = 0.4, label = T)
ovary <- RunUMAP(ovary, dims = 1:7)  #1:7/9/17
DimPlot(ovary, reduction = "umap", pt.size = 0.2, label = T)


# 筛选Marker
{
  VlnPlot(object = ovary, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) 
  # 找Oocyte的 cellmarker
  cluster3.markers <- FindMarkers(ovary, ident.1 = 3, ident.2 = c(0, 2,5), min.pct = 0.25)
  head(cluster3.markers, n = 60)
  cluster31.markers <- FindMarkers(ovary, ident.1 = 3, ident.2 = c(0,1,2,4,5,6), min.pct = 0.25)
  head(cluster31.markers, n = 60)
  
  VlnPlot(ovary, features = c("DCDC2","MYOCOS","NLRP14","NLRP5","RBM46","ENSMFAG00000032710",
                              "ZAR1","ZP2","PADI6","ENSMFAG00000046178","CLDN10","ZAR1L"))
  VlnPlot(ovary, features = c("RSPO2","BRDT","ESRP1","NLRP11","NLRP7","NLRP2","NLRP4",
                              "WEE2","DPPA3","TDRD1","NLRP13","LHX8"))  # Oo与其他cluster比较
  
  VlnPlot(ovary, features = c("BCL2L10","NPM2","ENSMFAG00000001960","KHDC3L","DPPA5","ACE2","AVEN",
                              "ESCO2","FBXW12","ELAVL2","LMOD3","KPNA7","H1FOO","OOEP","LDHAL6A","RNF17")) 
  VlnPlot(ovary, features = c("ENSMFAG00000043380","SLC25A31","ENSMFAG00000031694","IGF2BP3","TGM7","HENMT1",
                              "DLGAP5","PPP1R3A","AK5","FIGLA","CENPE","UHRF1","NEFM","SYCP2","MAEL","OTX2",
                              "ZP4","FKBP6","STK31","GTSF1"))  # Oo与其他cluster比较
  
  VlnPlot(ovary, features = c( "ACE2","ACTL8","ADAD1","ADAD2","AK5",
                               "ALPL","FHOD3","CENPO")) 
  
  # OO 5
  DotPlot(object = ovary, features = c("SYCP3","DDX4","GDF9","ZP3"))
  # GC 2 4
  DotPlot(object = ovary, features = c("AMH","NR5A2","CYP19A1","INHA","WT1"))
  # SC 0 11 13
  DotPlot(object = ovary, features = c("TCF21","STAR","ALDH1A1","COL1A1","COL1A2"))
  # SMC 1  14
  DotPlot(object = ovary, features = c("DES","RARB","CSRP1","ACTA2","ACTG2"))
  # EC 6
  DotPlot(object = ovary, features = c("CDH5","ERG","VWF","RNASE1"))
  # M 8 
  DotPlot(object = ovary, features = c("CD68","CD14","CD163"))
  # NKT 3
  DotPlot(object = ovary, features = c("CD3D","KLRB","REL","TRAC","CD3G","KLRD1"))  #KLRB报错
  
}
{
  #补充--用CellMarker验证
  DotPlot(ovary,feature=c("SYCP3","DDX4","GDF9","ZP3","FIGLA","LMOD3","RBM46","NETO1")) # 前5个为鉴定marker，后3个为新marker
  FeaturePlot(ovary, features = c("SYCP3","DDX4","GDF9","ZP3","FIGLA"))    # Oocyte       -- "SYCP3",  "DDX4","GDF9","ZP3",   "FIGLA"     @ Cluster 3
  FeaturePlot(ovary, features = c("LMOD3","RBM46","NETO1"))        # Oocyte novel --"LMOD3","RBM46","NETO1"
  
  VlnPlot(ovary, features = c("AMH","WT1","INHA","CYP19A1"))   # GC -- "AMH",  "WT1","INHA","CYP19A1",   "NR5A2"           @Cluster 1
  VlnPlot(ovary, features = c("TCF21","COL1A2","STAR","COL1A2","ALDH1A1"))  # Stromal Cell -- "TCF21",  "COL1A2","STAR"      @Cluster 0
  VlnPlot(ovary, features = c("DES","ACTA2","STAR","CSRP1","ACTG2"))    # Smooth Muscle Cells -- "DES",   "ACTA2","STAR",   "RARB"(有差别)
  VlnPlot(ovary, features = c("CDH5","VWF","ERG","RNASE1"))    # Endothelial cells-- "CDH5",   "VWF",  "ERG"            @ Cluster 4
  VlnPlot(ovary, features = c("CD3D","KLRB1","REL","CD3G","KLRD1"))   # Natural Killer T --"CD3D","KLRB1",  "REL"            @ Cluster 2
  VlnPlot(ovary, features = c("CD68","CD14","CD163"))     # Microphages-- "CD68",  "CD14"
  
}

# 筛选marker @240512
DimPlot(ovary, reduction = "pca") + NoLegend()
ovary.markers <- FindAllMarkers(ovary, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
# mk <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
rm(top10)
rm(mk)
top10 <- ovary.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 0.5) %>%  # 在每个分组中筛选出 avg_log2FC 大于5的行（avg_log2FC 大于5的结果已经很明显了）
  slice_head(n = 5) %>% # 从满足条件的结果中选取每个分组的前15条记录
  ungroup()
# features中有重复项，需删去
{unique_factors <- unique(top10$gene)
  print(unique_factors)
  top10 <- top10[!duplicated(top10$gene), ]
  top10$gene <- factor(top10$gene, levels = unique(top10$gene))
  }
mk <- top10$gene

#DoHeatmap(scRNA,features = mk)
DotPlot(ovary,features = mk)+RotatedAxis()

# As for oocyte
mk_oo <- c('LMOD3','RBM46','NETO1',  # novel gene identified by this article
           "SYCP3","DDX4","GDF9","ZP3", # mark gene
           "NLRP14","NLRP5","MYOCOS","NPM2","ZAR1" # selected this time
)
DotPlot(ovary,features = mk_oo)+RotatedAxis()  # 全都OK，特异性很好
DoHeatmap(ovary,features = mk_oo)
# ZAR1:zygote arrest 1--- This maternal effect gene is oocyte-specific and encodes a protein that is thought to function in the initiation of embryogenesis. 
#  A similar protein in mouse is required for female fertility

# 原文 Marker
rm(mk)
mk <- c("TCF21","COL1A2","STAR",  # SC
        "AMH","WT1","INHA","CYP19A1", # GC
        "DES","ACTA2",  # SMC
        "SYCP3","DDX4","GDF9","ZP3",# oocyte
        "CD3D","KLRB1",  # NKT
        "CDH5","VWF",  # EC
        "CD68","CD14")  # Macrophage
DotPlot(ovary,features = mk)+RotatedAxis()
DoHeatmap(ovary,features = mk)

new.cluster.ids <- c("SC","GC","SMC","Oocyte","NKT","EC","Macrophage")
names(new.cluster.ids) <- levels(ovary)
ovary <- RenameIdents(ovary, new.cluster.ids)
DimPlot(ovary, reduction = "tsne", label = TRUE, pt.size = 0.35) + NoLegend()
DimPlot(ovary, reduction = "umap", label = TRUE, pt.size = 0.35) + NoLegend()

# 在 ovary 的meta,data里加上 label   ovary@active.ident
ovary@meta.data$labels <- ovary@active.ident  # SC、GC、OO等注释信息

# 使用ifelse函数根据ovary@meta.data$orig.ident的不同分组为ovary@meta.data$stage添加了old和young的注释。
# 如果orig.ident的值在"group1"或"group3"中，那么age_group将被标记为"old"，否则为"young"。
ovary@meta.data$stage <- ifelse(ovary@meta.data$orig.ident %in% 
                                  c('OF1', 'OF2', 'OF3', 'OF4'), "old", 
                                "young")
# ovary@meta.data$stage <- NULL  


# 筛选后的oocyte数量还是少了点  @240512
saveRDS(ovary,file = 'ovary_tsne-240513-OK')

# ovary 分组: old+young OK --------
rm(list=ls())
gc()
setwd("/home/dengxy/scRNA-seq/代码复现/GSE130664-Cynomolgus Monkey Ovary/")
options(stringsAsFactors = F)
library(patchwork)
library(Seurat)
library(tidyverse)
library(data.table)
library(harmony)
library(pryr)  #install.packages('pryr')
mem_used()
# 读取细胞注释后的RDS文件
ovary <- readRDS("ovary_tsne-240513-OK")

ovary_young <- subset(ovary, subset = stage == 'young')
ovary_old <- subset(ovary, subset = stage == 'old')

ovary_young <- NormalizeData(ovary_young) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA()
ElbowPlot(ovary_young)

ovary_old <- NormalizeData(ovary_old) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>%
  RunPCA() 
ElbowPlot(ovary_old)

# Harmony 整合
ovary_young <- RunHarmony(ovary_young, group.by.vars = 'orig.ident')  # group.by.vars 参数提供的是元数据列的名称，而不是列中的具体值
ovary_old <- RunHarmony(ovary_old, group.by.vars = 'orig.ident') 

head(ovary_old@meta.data) # 查看元数据
mean(ovary_old$percent.mito) # Y-0.1961714  O-0.193435


DimPlot(ovary_young, reduction = 'pca')
DimPlot(ovary_old, reduction = 'pca')

ovary_old@meta.data$stage <- "old"
ovary_young@meta.data$stage <- "young"
class(ovary_young@meta.data$stage) # 查看加上的 label 是哪种数据类型-character

saveRDS(ovary_young,'ovary_young_harmony_240518.rds')
saveRDS(ovary_old,'ovary_old_harmony_240518.rds')

ovary_young <- readRDS("ovary_young_harmony_240518.rds") # 数据已整合
ovary_old <- readRDS('ovary_old_harmony_240518.rds')  # # 数据已整合

# ovary_young / old 细胞通讯 OK --------------------
{
  library(NMF)
  library(circlize)
  library(ComplexHeatmap)
  library(BiocNeighbors)
  library(dplyr)
  library(igraph)
  library(CellChat)
  library(patchwork)
  library(Seurat)
  options(stringsAsFactors = FALSE)
  
  ovary_young <- readRDS("ovary_young_harmony_240518.rds") # 数据已整合
  ovary_old <- readRDS('ovary_old_harmony_240518.rds')  # # 数据已整合
  
  data.input = LayerData(ovary_old, assay="RNA", layer='data')  # normalized data matrix
  meta = ovary_old@meta.data # a dataframe with rownames containing cell mata data
  cell.use = rownames(meta) # 检查一下是不是包括了 OF1,2,3,4和YF1,2,3,4
  
  data.input = data.input[, cell.use]  
  meta = meta[cell.use, ]  
  unique(meta$labels)
  
  cellchat <- createCellChat(object = data.input,meta = meta, group.by = "labels") # 这里的labels是此前细胞注释所定义的细胞类群，如SC GC SMC Oocyte等
  
  cellchat <- addMeta(cellchat, meta = meta)
  cellchat <- setIdent(cellchat, ident.use = "labels") 
  levels(cellchat@idents) 
  # "SC"         "GC"         "SMC"        "Oocyte"     "NKT"        "EC"         "Macrophage"
  groupSize <- as.numeric(table(cellchat@idents)) 
  
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  #考虑到食蟹猴与人类之间的亲缘关系较小鼠更近，因此这里采用human的数据进行分析
  showDatabaseCategory(CellChatDB)
  
  dplyr::glimpse(CellChatDB$interaction)  # Show the structure of the database
  CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat,features = NULL) #--- This step is necessary even if using the whole database
  future::plan("multisession", workers = 4)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)  # 略久
  cellchat <- projectData(cellchat, PPI.human) #投影到PPI（RNA到蛋白质）,运行时间略久
  # 计算通信概率并推断cellchat网络
  cellchat <- computeCommunProb(cellchat, raw.use = F)   # 运行时间超级久
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  
  cellchat <- netAnalysis_computeCentrality(cellchat, 
                                            slot.name = "netP")
  # 各细胞亚群之间的所有交互数量和强度
  groupSize <- as.numeric(table(cellchat@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(cellchat@net$count, 
                   vertex.weight = groupSize,  
                   weight.scale = T,   
                   label.edge= F,   
                   title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, 
                   vertex.weight = groupSize, 
                   weight.scale = T, 
                   label.edge= F, 
                   title.name = "Interaction weights/strength")
  
  # 各类细胞亚型所涉及的通路
  mat <- cellchat@net$weight
  par(mfrow = c(2,4), xpd=TRUE)  
  # mfrow 参数用于设置图形布局，指定将图形分成多少行和多少列
  # xpd 参数用于确定绘图是否可以扩展到设备的边界之外。当 xpd 设置为 TRUE 时，图形可以绘制到设备的边界之外
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), 
                   dimnames = dimnames(mat))  # 创建一个与原始矩阵 mat 相同大小的全零矩阵
    mat2[i, ] <- mat[i, ]  # 将矩阵 mat 中的第 i 行数据复制到新创建的矩阵 mat2 中
    netVisual_circle(mat2, vertex.weight = groupSize, 
                     weight.scale = T, 
                     edge.weight.max = max(mat), 
                     title.name = rownames(mat)[i])
    # 参数 vertex.weight、weight.scale、edge.weight.max 和 title.name 分别用于设置节点权重、权重缩放、边权重上限和图的标题。
  }
  cellchat@netP[["pathways"]]  # 查看一共有哪些通路
  # saveRDS(cellchat,'ovary_young_cellchat_240518_OK.rds')
  # saveRDS(cellchat,'ovary_old_cellchat_240518_OK.rds')