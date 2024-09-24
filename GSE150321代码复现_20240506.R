# 240501~240502 文章复现 @DXY -------

# Tissue 1 细胞分群注释 --------
rm(list=ls())   #清空环境变量
gc()

getOption('timeout')
options(timeout=10000000)
future::plan("multisession", workers = 4)

setwd("/home/dengxy/scRNA-seq/DAY2+DAY3/GSE150321_RAW/")
# Article: Cellular heterogeneity landscape in laryngeal squamous cell carcinoma

options(stringsAsFactors = F)
library(patchwork)
library(Seurat)
library(tidyverse)
library(data.table)
#library(harmony)
memory.limit()  #查看分配的内存大小
memory.size()  #所占内存
memory.limit(size = 9999990000000) #扩大内存

start = Sys.time()
# pd1里就有12 985个细胞,原文献也说只有12 985个细胞,
# 所以第二个csv不是另一个样本,不用进行读取

### 解压至csv文件，类似于matrix稀疏矩阵
pd1 <- fread('GSM4546857_LSCC01_DBEC_UMI.csv',sep=',',header = T,data.table=F) %>% column_to_rownames('Cell_Index')
# pd1 <- read.csv() #也可用read.csv()进行读取，速度较慢

# head(pd1[, 1:5])  # 原csv表格中，行名为细胞，列名为基因。后续转为Seurat对象时需要进行转置
# class(pd1) # 查看对象的类型  data.frame
# nrow(pd1)  # 获取行数  12985
# ncol(pd1)  # 获取列数  44808

pd1 <- as.data.frame(t(pd1))  #表达矩阵的读入格式：行名为基因，列名为细胞
head(pd1[, 1:5])
end = Sys.time()
end - start

# pd2 <- fread('GSM4546858_LSCC02_DBEC_UMI.csv',sep=',',header = T,data.table=F) %>% column_to_rownames('V1')
# pd2 <- as.data.frame(t(pd2))

# intergene <- intersect(rownames(pd1),rownames(pd2))
# pd1 <- pd1[intergene,]
# pd2 <- pd2[intergene,]

scRNA <- CreateSeuratObject(
  counts = cbind(pd1)) %>% 
  Seurat::NormalizeData() %>%  # 去除样本/细胞效应
  FindVariableFeatures(selection.method = "vst", nfeatures = 2500) %>% 
  ScaleData()  # 去除基因效应,在绘制热图的时候会需要使用
scRNA <- RunPCA(scRNA, npcs = 50, verbose = T) #npcs默认50
ElbowPlot(scRNA)
# PC控制在15以内，前10个PC已趋于缓和
# scRNA@meta.data$patient <- c(rep('LSCC1',ncol(pd1)),rep('LSCC2',ncol(pd2)))

# scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
print(scRNA[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scRNA, dims = 1:10, reduction = "pca")
DimPlot(scRNA, reduction = "pca")
DimHeatmap(scRNA, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(scRNA, dims = 1:15, cells = 500, balanced = TRUE)


#---质控
# Add number of genes per UMI for each cell to metadata
scRNA$log10GenesPerUMI <- log10(scRNA$nFeature_RNA) / log10(scRNA$nCount_RNA)  # 计算每个细胞的基因检测饱和度
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-") 
VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


# 滤过
mean(scRNA$nCount_RNA)  # 计算nCount_RNA的平均数
# 3077.34
median(scRNA$nCount_RNA)  # 计算nCount_RNA中位数
# 1946
table(scRNA$percent.mt<20)
table(scRNA$nFeature_RNA>200)  # 11138
table(scRNA$nCount_RNA>200) # 12985
table(scRNA$log10GenesPerUMI > 0.75) # 10862 
table(scRNA$nFeature_RNA>200 & scRNA$log10GenesPerUMI > 0.75 ) # 10637

table(scRNA$log10GenesPerUMI > 0.73 & scRNA$nFeature_RNA>200 & scRNA$nFeature_RNA<5000 & scRNA$nCount_RNA<22000) # 10700

# 原文筛选后保留了 10699 个细胞
# 筛选条件为  Cells were removed when the number of detected genes was less than 200

scRNA <- subset(scRNA, subset = log10GenesPerUMI > 0.73 &
                  nFeature_RNA>200  &
                  nFeature_RNA <5000 &
                  nCount_RNA<22000)     # 筛选到10700个细胞, 

scRNA <- FindNeighbors(scRNA,dims = 1:11) 
scRNA <- FindClusters(scRNA, resolution = 0.02) # 0.012~0.024
# 原文分了5个亚群
#   0     1    2    3    4  
# 5777  3873  493  366  190

table(scRNA$seurat_clusters)  #查看各聚类的细胞数量
# 0    1    2    3    4     # nFeature + UMI
# 5694 3869  511  368  195
# 0    1    2    3    4     # nFeature 
# 6226 3836  515  366  195 
prop.table(table(scRNA$seurat_clusters)) #查看各聚类的细胞比例
round(prop.table(table(scRNA$seurat_clusters)) * 100, 3)


# tSNE不影响各聚类的细胞比例和细胞数量   240504 @ DXY
scRNA <- RunTSNE(scRNA,dims = 1:20) #这一步很慢
# tSNE的 dims默认1:5, dims数量会影响聚类之间的分离效果，但对聚类后的细胞数量和比例无影响
DimPlot(scRNA,reduction = 'tsne')



scRNA.markers <- FindAllMarkers(scRNA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
# mk <- scRNA.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

rm(top10)
rm(mk)
top10 <- scRNA.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%  # 在每个分组中筛选出 avg_log2FC 大于5的行（avg_log2FC 大于5的结果已经很明显了）
  slice_head(n = 5) %>% # 从满足条件的结果中选取每个分组的前15条记录
  ungroup()
# features中有重复项，需删去
{unique_factors <- unique(top10$gene)
  print(unique_factors)
  top10 <- top10[!duplicated(top10$gene), ]
  top10$gene <- factor(top10$gene, levels = unique(top10$gene))
  str(scRNA)
  str(top10)
  }
mk <- top10$gene

#DoHeatmap(scRNA,features = mk)
DotPlot(scRNA,features = mk)+RotatedAxis()

# === === === === === ===
## 根据上述条件自己筛选的 marker及后续注释的细胞类型+相应基因
# "KRT19"  "KRT17"  "ABCC5"  "GSTP1"  "IGFBP2"   # Identity 0 Basal cell :KRT19,KRT17,KRT5
# "PTPRC"  "SRGN"   "IL32"   "CREM"   "CD2"      # Identity 1 Immune cell :均可，指向多种免疫细胞
# "WFDC2"  "CLU"    "SLPI"   "PIGR"   "KRT7"     # Identity 2 Epithelial cell :KRT7
# "SPARC"  "COL3A1" "COL1A2" "BGN"    "COL6A2"   # Identity 3 Fibroblast cell :COL3A1,COL1A2,COL6A2
# "AQP1"   "RNASE1" "PLVAP"  "VWF"    "EGFL7"    # Identity 4 Endothelial cell :AQP1,PLVAP,VWF,EGFL7
# === === === === === ===

# 原文给出的 marker
DotPlot(scRNA,features = c('KRT5',  # Tumor cells
                           'PTPRC', # Immune cells  筛到了
                           'CLU',   # Epithelial cells  筛到了 ，用 KRT7 更好
                           'AQP1',  # Endothelial cells  筛到了
                           'COL3A1','COL1A2')) +RotatedAxis()  # Fibroblasts  筛到了
FeaturePlot(scRNA, features = c('KRT5',  # Tumor cells
                                'PTPRC', # Immune cells  
                                'CLU',   # Epithelial cells  
                                'AQP1',  # Endothelial cells  
                                'COL3A1','CDKN2A'))  # Fibroblasts

# 把筛选到的marker和论文中给出的参考marker放一起比较,发现两种marker都可用
# 自己根据 avg_log2FC 筛选出来的marker ，平均表达量不如原文的 marker
# 在筛选指标里添加 group_by(cluster) 之后，筛出的marker多数与原文 marker一致，在细胞中的表达比例和平均表达量也符合筛选要求。  Very good!
DotPlot(scRNA,features = c("KRT19",  "KRT17",  "ABCC5",  "GSTP1",  "IGFBP2",  'KRT5', # Tumor cells  0
                           "SRGN",   "IL32",   "CREM",   "CD2",'PTPRC', # Immune cells  1
                           "WFDC2",   "SLPI",   "PIGR",   "KRT7" ,'CLU',  # Epithelial cells  2
                           "COL3A1", "COL1A2", "SPARC", "BGN", "COL6A2",  # Fibroblasts 3
                           "AQP1", "RNASE1", "PLVAP", "VWF", "EGFL7" # Endothelial cells 4
)) +RotatedAxis() 
# 细胞类型注释
new.cluster.ids <- c("Tumor cells", "Immune cells", "Epithelial cells","Fibroblasts cells", "Endothelial cells") # 因为自己分析聚类后得到的细胞类型与原文相似，仅identity 0 怀疑是Basal cell。故此处仍沿用文献中的注释名
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
scRNA@meta.data$labels <- scRNA@active.ident
DimPlot(scRNA, reduction = "tsne", label = TRUE, pt.size = 0.1) + NoLegend() 


# 画饼图
# 获取聚类结果和对应的细胞数量
cluster_table <- table(Idents(scRNA))
# 创建一个数据框，包含各聚类的名称和细胞数量
cluster_df <- data.frame(Cluster = names(cluster_table), Count = as.vector(cluster_table))
# 计算各聚类的百分比
cluster_df$Percentage <- cluster_df$Count / sum(cluster_df$Count) * 100
# 绘制饼状图和添加百分比标签
library(ggplot2)
ggplot(cluster_df, aes(x = "", y = Percentage, fill = Cluster)) +  # x 轴为空（因为饼图没有实际的 x 轴）
  geom_bar(stat = "identity", width = 1) +  # 添加柱状条，stat = "identity" 表示使用数据中的实际值作为条的高度
  coord_polar("y", start = 0) +  # 将坐标系转换为极坐标系，使得柱状条形成一个饼图。"y" 参数表示在 y 轴上进行极坐标转换，start = 0 表示从 0 度开始绘制
  geom_text(aes(label = paste0(round(Percentage), "%")), # 将百分比数据转换为字符型，并添加 % 符号
            position = position_stack(vjust = 0.5)) +  # 指定了文本标签的位置，vjust = 0.5表示文本在垂直方向上居中对齐
  labs(fill = "Cluster") +  # 更改图例的标题，将其改为 “Cluster”
  theme_void() +  # 移除默认的背景、网格线等元素，使得画布是空白的
  theme(legend.position = "right")  # 设置图例的位置

saveRDS(scRNA,file = 'LSCC_tissue1_240507.RDS')


# Tissue 1- Immune cells 亚群分类 --------
rm(list=ls())   #清空环境变量
gc()
library(Seurat)
setwd("/home/dengxy/scRNA-seq/DAY2+DAY3/GSE150321_RAW/")
scRNA <- readRDS('LSCC_tissue1_240507.RDS')
Immune <- subset(scRNA, idents = 'Immune cells')

Immune <- NormalizeData(Immune) %>%
  FindVariableFeatures() %>%
  ScaleData() 

Immune <- RunPCA(Immune,npcs = 50, verbose = T)  
ElbowPlot(Immune) # 10PC时趋于平缓(不运行RunPCA)   # 10PC时趋近缓和（运行PCA）

Immune <- FindNeighbors(Immune,dims = 1:15) 
Immune <- FindClusters(Immune, resolution = 0.28)#resolution调分辨率

table(Immune$seurat_clusters)  #查看各聚类的细胞数量
prop.table(table(Immune$seurat_clusters)) #查看各聚类的细胞比例
round(prop.table(table(Immune$seurat_clusters)) * 100, 3)
# 原文分了8个亚群
# 1208    791    783     454   291  144  114  88

Immune.markers <- FindAllMarkers(Immune, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
rm(top10)
rm(mk)
top10 <- Immune.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%  # 在每个分组中筛选出 avg_log2FC 大于5的行（avg_log2FC 大于5的结果已经很明显了）
  slice_head(n = 5) %>% # 从满足条件的结果中选取每个分组的前15条记录
  ungroup()
# features中有重复项，需删去
{unique_factors <- unique(top10$gene)
  print(unique_factors)
  top10 <- top10[!duplicated(top10$gene), ]
  top10$gene <- factor(top10$gene, levels = unique(top10$gene))
  str(Immune)
  str(top10)
  }
mk <- top10$gene

# DoHeatmap(Immune,features = mk)

# 原文 marker
mk <- c('CD3D',          # T
        'IL7R','CXCR4',  #0 Native T             
        'FOXP3',         #1 Regulatoty T        
        'CXCL13',        #3 Helper T
        'GZMA','GZMB','GNLY',  #2 Cytotoxic T   
        'MKI67','HMGB2', #7 Progenitor T        
        'MS4A1',         #5 B                   
        'IGHG3','IGHG4', #6 Plasma cell(B)      
        'APOC1',         #4 Macrophage          
        'CD163')         #4 M2 macrophage       


DotPlot(Immune,features = mk) + RotatedAxis()
VlnPlot(Immune,features = mk)
#DoHeatmap(Immune,features = mk)
VlnPlot(Immune,features = c('CD3D', 'MS4A1','APOC1'))
VlnPlot(Immune,features = c('CD3D',  # T
                            'IL7R','CXCR4', # Native T              
                            'FOXP3',  # Regulatoty T               
                            'CXCL13', # Helper T
                            'GZMA','GZMB','GNLY',  # Cytotoxic T   
                            'MKI67','HMGB2'))  # Progenitor T
VlnPlot(Immune,features = c('CXCR4','LMNA','MYADM'))  # Native T 
VlnPlot(Immune,features = c('FOXP3','IL2RA','BATF'))  # Regulatoty T 

FeaturePlot(Immune, features = c("CD163", "CXCL8"))

Immune <- RunTSNE(Immune,dims = 1:20) 
DimPlot(Immune,reduction = 'tsne')

# 细胞类型注释
new.cluster.ids <- c("Cytotoxic T","Regulatoty T","Helper T",'Macrophage',
                     "Cytotoxic T-2",'B','Progenitor T','Plasma cell') 
# new.cluster.ids <- c('0','1','2','3','4','5','6','7')
names(new.cluster.ids) <- levels(Immune)
Immune <- RenameIdents(Immune, new.cluster.ids)
DimPlot(Immune, reduction = "tsne", label = TRUE, pt.size = 0.2) + NoLegend() 

# GO分析
library(IRanges)
library(S4Vectors)
library(AnnotationDbi)
library(stats4)
library(BiocGenerics)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(enrichplot)
rm(diff)
rm(sig_diff)

new.cluster.ids <- c('0','1','2','3','4','5','6','7')
names(new.cluster.ids) <- levels(Immune)
Immune <- RenameIdents(Immune, new.cluster.ids)
diff <- FindMarkers(Immune, ident.1 = 5, ident.2 = c(0,1,2,3,4))  #在Immune中，5是B cell
sig_diff <- diff %>% 
  filter(p_val_adj<0.05 & abs(avg_log2FC)>0.25) %>% 
  slice_head(n = 10)

#差异基因GO富集分析
ego_ALL <- enrichGO(gene          = row.names(sig_diff),
                    #universe     = row.names(sig_diff),
                    OrgDb         = 'org.Hs.eg.db',
                    keyType       = 'SYMBOL',
                    ont           = "ALL",
                    pAdjustMethod = "BH",
                    pvalueCutoff  = 0.01,
                    qvalueCutoff  = 0.05)
ego_all <- data.frame(ego_ALL)

dotplot(ego_ALL,split="ONTOLOGY")+
  facet_grid(ONTOLOGY~.,scales = "free")#可视化



# Tissue 1- Tumor cells 亚群分类 --------
Tumor <- subset(scRNA, idents = 'Tumor cells')
Tumor <- NormalizeData(Tumor) %>%
  FindVariableFeatures()  %>%
  ScaleData()  %>%
  RunPCA() 

ElbowPlot(Tumor) # 15~20 PC时趋于平缓（不运行RunPCA）   # 10PC时趋于平缓（运行RunPCA）

Tumor <- FindNeighbors(Tumor,dims = 1:20)
Tumor <- FindClusters(Tumor, resolution = 0.3)#resolution调分辨率
Tumor$RNA_snn_res.0.25 <- NULL

table(Tumor$seurat_clusters)  #查看各聚类的细胞数量
prop.table(table(Tumor$seurat_clusters)) #查看各聚类的细胞比例
# 原文分了5个亚群
# 3374  1008  948  383  64

Tumor.markers <- FindAllMarkers(Tumor, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
library(dplyr)
rm(top10)
rm(mk)
top10 <- Tumor.markers %>%
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


mk <- c('SPRR1A','SPRR1B',  # Keratinocyte cancer   4
        'BIRC5','TPX2',  # Proliferative cancer     3
        'SERPINB5','AMD1',  # Metastasis cancer     1
        'IL32','CD2',  # Immune cancer              5
        'MTRNR2L1','MTRNR2L8')  # Immortal cancer   0+2

DotPlot(Tumor,features = mk) + RotatedAxis()
VlnPlot(Tumor,features = mk)
#FeaturePlot(Tumor,features = mk)
DoHeatmap(Tumor,features = mk)

Tumor <- RunUMAP(Tumor,dims = 1:15)
DimPlot(Tumor,reduction = 'umap')
Tumor <- RunTSNE(Tumor,dims = 1:20) 
DimPlot(Tumor,reduction = 'tsne')

FeaturePlot(Tumor,features = c('SPRR3','SPRR2A','SPRR1A','SPRR2D'))  # Keratinocyte cancer
FeaturePlot(Tumor,features = c('MKI67','TPX2','CKS1B','CDCA3'))      # Proliferative cancer 



# 细胞类型注释
new.cluster.ids <- c( "Immortal cancer-1","Metastasis cancer", "Immortal cancer-2","Proliferative cancer", "Keratinocyte cancer",'Immune cancer') 
names(new.cluster.ids) <- levels(Tumor)
Tumor <- RenameIdents(Tumor, new.cluster.ids)
Tumor@meta.data$tsne_celltype <- Tumor@active.ident
DimPlot(Tumor, reduction = "tsne", label = TRUE, pt.size = 0.2) + NoLegend() 

table(Tumor$percent.mt >50)
VlnPlot(Tumor,features = 'percent.mt')  # Immortal cancer-1,2 的线粒体基因异常高，这里的线粒体基因没有过滤，特此注明。

saveRDS(Tumor,'Tumor_tsne_240606.rds')
