## scRNA-seq 数据处理 ----------
rm(list=ls())   #清空环境变量
gc()
future::plan("multisession", workers = 4)

setwd("/home/dengxy/scRNA-seq/代码复现/GSE255690_human_ovary/scRNA-seq/scRNA/")
options(stringsAsFactors = F)
library(patchwork)
library(Seurat)
library(tidyverse)
library(data.table)
library(harmony)

## 前期 ------------预处理文件--------------------
# 使用R语言中的file.rename函数来批量处理文件名，将目录中所有文件名中的“-”替换为“_”
directory <- "/home/dengxy/scRNA-seq/代码复现/GSE255690_human_ovary/scRNA-seq/scRNA/"  
files <- list.files('./','^GSM')  
for (file in files) {
  # 检查文件名中是否包含"-"
  if (grepl("-", file)) {
    # 新文件名，将"-"替换为"_"
    new_file <- gsub("-", "_", file)
    # 获取完整的文件路径
    old_filepath <- file.path(directory, file)
    new_filepath <- file.path(directory, new_file)
    # 重命名文件
    file.rename(old_filepath, new_filepath)
    cat("Renamed:", file, "->", new_file, "\n")
  }
}

fs=list.files('./','^GSM')  # 列出当前目录下所有开头是GSM的文件
#e.g. GSM8077839_st_old-5_matrix.mtx.gz

# 然后获取四个样本信息
library(stringr)
samples=str_split(fs,'_',simplify = T)[,1]  #以‘_’为分界线，保留隔开后的第一个字符串
# e.g. GSM8077839
lapply(unique(samples),function(x){
  y=fs[grepl(x,fs)]
  folder=paste(str_split(y[1],'_',simplify = T)[,1:3],collapse = '')  
  dir.create(folder,recursive = T)
  file.rename(y[1],file.path(folder,"barcodes.tsv.gz"))  # 细胞
  file.rename(y[2],file.path(folder,"features.tsv.gz"))  # 基因
  file.rename(y[3],file.path(folder,"matrix.mtx.gz"))    # 稀疏矩阵
})


## 前期 ----------  批量读取文件----------------------
# 设置路径循环
assays <- dir('/home/dengxy/scRNA-seq/代码复现/GSE255690_human_ovary/scRNA-seq/scRNA') 
dir <- paste0('/home/dengxy/scRNA-seq/代码复现/GSE255690_human_ovary/scRNA-seq/scRNA/',assays) #注意这行末尾有‘/’
# 批量读取数据
scRNAlist <- list()
# 如何获取分组信息？
sample_name <- c('young 1','young 3','young 4','middle 2','middle 3','middle 4','old 1','old 2','old 3')

for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts,project = sample_name[i],
                                       min.cells=3,min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]],add.cell.id=sample_name[i])
  if(T){
    scRNAlist[[i]][['percent.mt']] <- PercentageFeatureSet(scRNAlist[[i]], pattern = '^MT-')
  }
}


# 数据整合 
ovary <- merge(scRNAlist[[1]],scRNAlist[2:length(scRNAlist)])  # 26164x93217
ovary_young <- merge(scRNAlist[[1]],scRNAlist[2:3])
ovary_middle <- merge(scRNAlist[[4]],scRNAlist[5:6])
ovary_old <- merge(scRNAlist[[7]],scRNAlist[8:9])

setwd('/home/dengxy/scRNA-seq/代码复现/GSE255690_human_ovary/scRNA-seq/')
future::plan("multisession", workers = 4)


## 前期 ---各组ovary去批次---------------------------------------------
ovary_young <- NormalizeData(ovary_young) %>% 
  FindVariableFeatures() %>%
  ScaleData(ovary) %>%
  RunPCA(verbose = FALSE)  # verbose 参数控制是否在运行过程中输出详细信息

mk <- c('DCN','STAR',  # Theca and Stroma cells--T&S
        'ACTA2','MUSTN1', # Smooth muscle cells (SMCs)
        'TM4SF1','VWF', # Endothelial cells(EC)
        'CCL5','NKG7', # Natural killer cells (NK)
        'IL7R','KLRB1',  # T lymphocytes
        'TYROBP', 'IFI30', # Monocytes
        'GSTA1','AMH','HSD17B1', # Granulosa cells(GCs)
        'ZP3','FIGLA'    # Oocytes  'TUBB8',
)
DotPlot(ovary,features = mk)+RotatedAxis()


ovary_young <- RunHarmony(ovary_young, group.by.vars = 'orig.ident' )  # group.by.vars 参数提供的是元数据列的名称，而不是列中的具体值
head(ovary_young@meta.data) # 查看元数据
mean(ovary_young$percent.mt) # Y-9.845159  

# saveRDS(ovary_young,'ovary_young.rds')   #数据已整合

#new.cluster.ids <- c("T&S","SMC","EC","Immune cell","GC","Oocyte") # resolution = 0.02
new.cluster.ids <- c("T&S","EC","SMC","Immune cell","monocyte","GC",'Stroma GC',"Unknown GC","Oocyte") # resolution = 0.1
#new.cluster.ids <- c(0,1,2,3,4,5,6,7,8)
names(new.cluster.ids) <- levels(ovary)
ovary <- RenameIdents(ovary, new.cluster.ids)
ovary@meta.data$umap_celltype <- ovary@active.ident
ovary <-RunUMAP(ovary, reduction = "harmony", dims = 1:25)
ovary <-RunTSNE(ovary, reduction = "harmony", dims = 1:20)

DimPlot(ovary, reduction = "umap",
        group.by = "umap_celltype",  # 根据指定的元数据列进行分组
        label = T,  # 是否在图中显示分组标签
        repel = T,  # 是否使用力导向算法避免标签重叠
        seed = 10,   # 随机数种子，用于打乱细胞顺序时保证结果可重复。默认为1
        pt.size = 0.35) + NoLegend()