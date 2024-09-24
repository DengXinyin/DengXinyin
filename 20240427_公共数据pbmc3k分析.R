### 根据B细胞和周围其他细胞类型之间的差异基因，进一步筛选B细胞亚群之间的差异基因
pbmc <- readRDS(file = "/home/dengxy/scRNA-seq/DAY4/pbmc3k-cellchat/pbmc3k_final-240119.rds") # 此前已经分析好的 pbmc 数据
DimPlot(pbmc,reduction = 'umap',label = T,repel = T)

library(Seurat)
pbmc_B <- subset(pbmc,seurat_annotations %in% 'B')  # 特定聚类，如这里为 B cell
pbmc_B <- NormalizeData(pbmc_B) %>%
  ScaleData() %>%
  FindVariableFeatures() %>%
  RunPCA()
ElbowPlot(pbmc_B)

pbmc_B <- FindClusters(pbmc_B, resolution = 1) # 4 clusters

pbmc_B <- RunUMAP(pbmc_B, dims = 1:10)
DimPlot(pbmc_B, reduction = "umap",group.by = 'seurat_clusters')
pbmc_B <- RunTSNE(pbmc_B, dims = 1:20)
DimPlot(pbmc_B, reduction = "tsne",group.by = 'seurat_clusters')

B.cell.markers <- FindMarkers(pbmc,ident.1 = '3',group.by = 'seurat_clusters') # cluster之前注释为B cell
library(dplyr)
rm(top5)
top5 <- B.cell.markers %>% 
  dplyr::filter(avg_log2FC > 5) %>%  # 在每个分组中筛选出 avg_log2FC 大于5的行（avg_log2FC 大于5的结果已经很明显了）
  slice_head(n = 5) %>% # 从满足条件的结果中选取每个分组的前5条记录
  ungroup()
mk <- rownames(top5)
mk
mk <- c("CD79A","MS4A1","LINC00926","TCL1A","VPREB3")   

DotPlot(pbmc, features = mk)
DotPlot(pbmc_B, features = mk) + NoLegend() 