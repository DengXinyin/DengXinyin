# 创建Seurat对象-------------------------
library(Seurat)  # 用于处理单细胞RNA测序数据
library(stringr)  # 用于字符串操作
library(fs)  # 用于文件系统操作
library(doParallel)  # 用于并行计算
library(foreach)  # 用于并行循环
setwd('/home/hsinyinteng/GSE260685_GSE260686_human_ovary/GSE260685_scRNA_seq/')

# 获取当前目录中的所有文件夹
dirs <- list.dirs(recursive = FALSE)  # 仅获取一级子目录

# 注册并行后端
Miss_Qiu_596 <- makeCluster(detectCores() - 15)  # 创建并行集群，使用所有可用核心减一
registerDoParallel(Miss_Qiu_596)  # 注册并行后端

# 并行处理每个文件夹并返回结果
seurat_objects_list <- foreach(dir = dirs, .packages = c('Seurat', 'Matrix')) %dopar% {
  # 读取数据并创建Seurat对象
  data <- Read10X(data.dir = dir)  # 从目录中读取10X Genomics数据
  seurat_object_name <- basename(dir)  # 获取文件夹名作为Seurat对象的名称
  seurat_object <- CreateSeuratObject(counts = data, 
                                      assay = "RNA", 
                                      project = seurat_object_name) 
  mt.genes <- grep(pattern = "^MT-", x = rownames(seurat_object), value = TRUE)  
  rownames(seurat_object@assays$RNA@layers$counts) <- rownames(seurat_object)
  colnames(seurat_object@assays$RNA@layers$counts) <- colnames(seurat_object)
  seurat_object$percent.mito <- (Matrix::colSums(seurat_object@assays$RNA@layers$counts[mt.genes, ])/Matrix::colSums(seurat_object@assays$RNA@layers$counts))*100
  
  # 返回Seurat对象列表
  list(name = seurat_object_name, object = seurat_object)
  
}

# 停止并行集群
stopCluster(Miss_Qiu_596)  # 停止并行集群

# 将返回的Seurat对象分配到全局环境中，并重新命名对象
seurat_objects <- list()
for (obj in seurat_objects_list) {
  assign(obj$name, obj$object, envir = .GlobalEnv)
  seurat_objects[[obj$name]] <- obj$object
}

# 打印seurat_objects列表内容以验证
print(seurat_objects)

# 保存Seurat对象到文件
# saveRDS(seurat_objects, "GSE260685_scRNA_seq_20240719.rds")

# 输出创建的Seurat对象名称
print(ls(pattern = "^s"))  # 打印所有以“s”开头的对象名称
rm(seurat_objects_list)
rm(obj)

rm(Miss_Qiu_596)

# 结合文献进行质控  ------------------
setwd('/home/hsinyinteng/GSE260685_GSE260686_human_ovary/GSE260685_scRNA_seq/')

library(Seurat)
GSE260685_spatial_objects_all <- readRDS('GSE260685_scRNA_seq_20240719.rds')

# # 如果前面没有计算线粒体基因
# mt.genes <- grep(pattern = "^MT-", x = rownames(s3), value = TRUE)  
# rownames(s3@assays$RNA@layers$counts) <- rownames(s3)
# colnames(s3@assays$RNA@layers$counts) <- colnames(s3)
# 
# s3$percent.mito <- (Matrix::colSums(s3@assays$RNA@layers$counts[mt.genes, ])/Matrix::colSums(s3@assays$RNA@layers$counts))*100
# # genes_to_keep <- setdiff(names(which(Matrix::rowSums(s3@assays$RNA@layers$counts )>5)),mt.genes)
# # s3_subset <- subset(s3, features =genes_to_keep, subset = nFeature_Spatial > 300 & percent.mito < 30)
# 
# # 提取 nFeature_RNA 和 nCount_RNA 列并计算平均值
# metadata <- s3@meta.data  # 提取 Seurat 对象的元数据
# mean(metadata$nFeature_RNA, na.rm = TRUE)  # 1133.269
# mean(metadata$nCount_RNA, na.rm = TRUE)    # 3178.423
# mean(s3$percent.mito, na.rm = TRUE)        # 13.56081
# table(s3$nCount_RNA > 2300 & s3$percent.mito < 15) 
# # FALSE  TRUE 
# # 6733  7600

s3 <- GSE260685_spatial_objects_all$s3
s4 <- GSE260685_spatial_objects_all$s4
s5 <- GSE260685_spatial_objects_all$s5

# 质控筛选细胞
s3_subset <- subset(s3, subset = nCount_RNA > 2300 & percent.mito < 15)  #36601x7600, 原文 7571 cells
s4_subset <- subset(s4, subset = nCount_RNA > 1275 & percent.mito < 12)  #36601x7340, 原文 7228 cells
s5_subset <- subset(s5, subset = nCount_RNA > 1096 & percent.mito < 15)  #36601x6361, 原文 6339 cells

# 原文描述如下：
# 来自三个供体的ScRNA数据最初由密歇根大学高级基因组学核心使用Cell Ranger v4.0.0进行处理。
# 主要步骤包括从原始双端测序读数中提取细胞条形码和 UMI、与人类 Ensembl 基因比对以及基于 UMI 的重复数据删除，从而产生逐个细胞的 UMI 计数表，
# 分别由三个样本的 14,322、13,901 和 9149 个细胞的“filtered_feature_bc_matrix”表示，共 20,886 个基因。
# 细胞过滤使用 （i） 最少数量的 UMI（在 Seurat 中称为“nCount”）和 （ii） 对应于线粒体编码基因的转录本百分比（“% MT”）。
# 三个样本的临界值各不相同：nCount：2300、1275 和 1096;% MT：分别为 15、12 和 15。
# 临界值是根据每个样本的 nCount 和 % MT 分布（未显示）选择的。
# 过滤后，还剩下 7571、7228 和 6339 个细胞供进一步分析。
# 这些细胞的 nCount 平均值、检测到的基因数 （“nFeature”） 和 % MT 如图 S1A 所示。

# s3+s4+s5 样本整合+去批次+细胞注释 --------------------
# 加载所需包
library(Seurat)
library(harmony)
future::multisession(4)

# 标记样本来源
# s3_subset$sample <- "s3"
# s4_subset$sample <- "s4"
# s5_subset$sample <- "s5"

# 合并 Seurat 对象
GSE260685_scRNA_all <- merge(s3_subset, 
                             y = list(s4_subset, s5_subset), 
                             add.cell.ids = c("s3", "s4", "s5"))
# 标准化和识别变异特征基因
GSE260685_scRNA_all <- NormalizeData(GSE260685_scRNA_all)  %>%
  FindVariableFeatures()  %>%
  ScaleData()  %>%
  RunPCA(npcs = 30)

# 使用 Harmony 去除批次效应
GSE260685_scRNA_all <- RunHarmony(GSE260685_scRNA_all, 
                                  "orig.ident", 
                                  plot_convergence = TRUE,
                                  reduction.use = "pca",
                                  dims.use = NULL,
                                  reduction.save = "harmony")



ElbowPlot(GSE260685_scRNA_all)  # 前15 PCs

GSE260685_scRNA_all <- FindNeighbors(GSE260685_scRNA_all, reduction = "harmony", dims = 1:15)
GSE260685_scRNA_all <- FindClusters(GSE260685_scRNA_all, resolution = 0.06) # 需要分成 4 clusters
GSE260685_scRNA_all$RNA_snn_res.0.045 <- NULL
#1:20+0.01
#1:15+0.05  /0.06***
#1:13+0.06* /0.05/0.04/0.03/0.02  /0.01
#1:12+0.01  /0.02  /0.03*
#1:11+0.01  /0.02  /0.03/0.04/0.05/0.06*
#1:10+0.01  /0.02  /0.03/0.04/0.05/0.06**

GSE260685_scRNA_all <- RunTSNE(GSE260685_scRNA_all, reduction = "harmony", dims = 1:30)
DimPlot(GSE260685_scRNA_all, reduction = "tsne", group.by = "orig.ident")
DimPlot(GSE260685_scRNA_all, reduction = "tsne", group.by = "seurat_clusters")

# 原文marker
mk <- c('DCN','PDGFRA','APOE','LUM','ARX','FHL2',  # Stromal cells-基质细胞
        'RGS5','NOTCH3','ACTA2','MUSTN1','PDGFRB',  # Pericytes-周细胞
        'VWF','PECAM1','CD34','CLDN5','NOTCH4',  # Endothelial cells-内皮细胞
        
        # 4 immune cell subtypes-免疫细胞亚型
        'CD3D','CD3G','IL7R',  # T cells
        'KLRD1','TRCD','GNLY',  # Natural Killer cells (NK)-自然杀伤细胞
        'CD68','CD14','FOLR2',  # Macrophages-巨噬细胞
        'KIT','TPSB2','TPSAB1'  # Mast cells-肥大细胞
)
DotPlot(GSE260685_scRNA_all,features = mk)+RotatedAxis()

GSE260685_scRNA_all <- RunUMAP(GSE260685_scRNA_all, reduction = "harmony", dims = 1:30)
DimPlot(GSE260685_scRNA_all, reduction = "umap", group.by = "orig.ident")

# 可视化批次效应移除的效果
# 整合前的 UMAP 图
# GSE260685_scRNA_all_pre_harmony <- ScaleData(GSE260685_scRNA_all, verbose = FALSE)
# GSE260685_scRNA_all_pre_harmony <- RunPCA(GSE260685_scRNA_all_pre_harmony, npcs = 30, verbose = FALSE)
# GSE260685_scRNA_all_pre_harmony <- RunUMAP(GSE260685_scRNA_all_pre_harmony, reduction = "pca", dims = 1:30, verbose = FALSE)
# DimPlot(GSE260685_scRNA_all_pre_harmony, reduction = "umap", group.by = "sample") + ggtitle("UMAP before Harmony")


# 连接数据层
GSE260685_scRNA_all <- JoinLayers(GSE260685_scRNA_all)

# 识别每个聚类的标志基因
markers <- FindAllMarkers(GSE260685_scRNA_all, 
                          only.pos = TRUE, 
                          min.pct = 0.25, 
                          logfc.threshold = 0.25)

# 查看标志基因
head(markers)

# 手动注释聚类
# 需要根据已知的细胞类型标志基因列表对每个聚类进行注释

# 为每个聚类添加注释
new_cluster_ids <- c("Stromal cells",
                     "T + NK",
                     "Macrophages",
                     "Endothelial + Mast cells", 
                     "Pericytes")
names(new_cluster_ids) <- levels(GSE260685_scRNA_all)

# 将注释添加到对象中
GSE260685_scRNA_all <- RenameIdents(GSE260685_scRNA_all, new_cluster_ids)
GSE260685_scRNA_all$celltype <- GSE260685_scRNA_all@active.ident
# 绘制带注释的 UMAP 图
DimPlot(GSE260685_scRNA_all, reduction = "tsne", label = TRUE, label.size = 5, group.by = "ident")


# 显示每个聚类的标志基因热图
library(dplyr)
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(GSE260685_scRNA_all, 
          features = top10$gene,
          group.by = "ident",
          group.bar = TRUE,
          group.colors = NULL,
          label = TRUE,
          size = 5.5,
          hjust = 0,
          vjust = 0,
          angle = 45,
          raster = TRUE,
          draw.lines = TRUE,
          lines.width = NULL,
          group.bar.height = 0.02,
          combine = TRUE) + NoLegend()

saveRDS(GSE260685_scRNA_all, 'GSE260685_scRNA_tSNE_20240721.rds')