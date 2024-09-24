## 整体大亚群注释 OK ------------
ovary <- readRDS('ovary_harmony_filted_240531.rds')
ovary <- FindNeighbors(ovary, reduction = "harmony", dims = 1:20) 
ovary <- FindClusters(ovary, resolution = 0.02)  # dims = 1:20, Y+M+O:0.01,0.02, 0.1  # Y-0.01  M-0.07  O-0.14
# 先用GC和Oocyte的marker看一下这两类细胞分布在哪些亚群
# GC和Oocyte marker来源：
# 结合人类卵巢marker@王世宣老师 Wu Meng et al., Nature Aging. 2024 ）+
#     食蟹猴卵巢marker@刘光慧老师（ Wang Si et al., Cell. 2020 ）+
#     人类卵巢marker@乔杰老师（Zhang Yaoyao et al.,Molecular cell. 2018）
# 注意：食蟹猴卵巢marker不一定适用于人类，仅作参考
DotPlot(ovary,features = c('GSTA1','AMH','HSD17B1', 
                           "INHA",'CYP11A1','STAR','TST', # granulosa cells(GCs)      
                           
                           'TUBB8','ZP3','FIGLA', 
                           "SYCP3","DDX4","GDF9",
                           'ZP2','SOX30','ZAR1',
                           'DAZL','YBX2','LHX8','NOL4'  # oocytes
))+RotatedAxis()
# 'LCP2'-oocyte
# 'ZEB2'-GC

# 人类卵巢marker@王世宣
rm(mk)  
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

# 结合人类卵巢marker@王世宣老师 Wu Meng et al., Nature Aging. 2024 ）+
#     食蟹猴卵巢marker@刘光慧老师（ Wang Si et al., Cell. 2020 ）+
#     人类卵巢marker@乔杰老师（Zhang Yaoyao et al.,Molecular cell. 2018）
# 注意：食蟹猴卵巢marker不一定适用于人类，仅作参考
rm(mk)
mk <- c('DCN','STAR',"TCF21","COL1A2",  # theca and stroma cells (T&S) 
        'ACTA2','MUSTN1',"DES", # smooth muscle cells (SMCs)
        'TM4SF1','VWF',"CDH5", # endothelial cells (EC)  
        
        # Immune
        'CCL5','NKG7', # natural killer cells (NK) 
        "CD3D",'IL7R','KLRB1','LAG3', 
        'TYROBP', 'IFI30', # monocytes
        "CD68","CD14",'CD4', # Macrophages
        
        'GSTA1','AMH','HSD17B1', # granulosa cells (GCs) 
        "INHA",'CYP11A1','TST',     # T&S-GC  'STAR',
        
        'TUBB8','ZP3','FIGLA',   # oocytes
        "SYCP3","DDX4","GDF9",
        'ZP2','SOX30','ZAR1','DAZL','YBX2','LHX8',
        'NOL4','GPD1','NTF4'
)
DotPlot(ovary,features = mk)+RotatedAxis()
DimPlot(ovary, reduction = "umap",label = T, pt.size = 0.35)
ovary@meta.data$RNA_snn_res.0.1 <- NULL

new.cluster.ids <- c("T&S","SMC","EC","Immune cell","GC","Oocyte") # resolution = 0.02
#new.cluster.ids <- c("T&S-1","SMC","EC","T&S-2","Immune cell","monocyte","GC",'T&S-GC',"Oocyte") # resolution = 0.1
#new.cluster.ids <- c(0,1,2,3,4,5,6,7,8)
names(new.cluster.ids) <- levels(ovary)
ovary <- RenameIdents(ovary, new.cluster.ids)
ovary@meta.data$umap_celltype <- ovary@active.ident  # SC、GC、OO等注释信息

ovary <-RunUMAP(ovary, reduction = "harmony", dims = 1:25)
ovary <-RunTSNE(ovary, reduction = "harmony", dims = 1:20)

DimPlot(ovary, reduction = "umap",
        group.by = "umap_celltype",  # 根据指定的元数据列进行分组
        label = T,  # 是否在图中显示分组标签
        repel = T,  # 是否使用力导向算法避免标签重叠
        seed = 10,   # 随机数种子，用于打乱细胞顺序时保证结果可重复。默认为1
        pt.size = 0.35) + NoLegend()

saveRDS(ovary,file='/home/dengxy/scRNA-seq/代码复现/GSE255690_human_ovary/ovary_harmony_filted_240531.rds')  # OK