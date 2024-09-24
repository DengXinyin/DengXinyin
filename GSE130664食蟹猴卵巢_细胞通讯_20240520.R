# ovary_young + ovary_old 细胞通讯 ---OK-----------------

  # Comparison analysis of multiple datasets using CellChat
  # Load the required libraries------
  rm(list=ls())
  gc()
  setwd("/home/dengxy/scRNA-seq/代码复现/GSE130664-Cynomolgus Monkey Ovary/")
  
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
  future::plan("multisession", workers = 4)
  
  # Load CellChat object of each dataset and then merge together-------
  # saveRDS(cellchat,'ovary_young_cellchat_240518_OK.rds')   #数据已完成细胞通讯分析
  # saveRDS(cellchat,'ovary_old_cellchat_240518_OK.rds')     #数据已完成细胞通讯分析
  ovary_young <- readRDS("ovary_young_cellchat_240518_OK.rds") 
  ovary_old <- readRDS('ovary_old_cellchat_240518_OK.rds')  
  ovary_young <- updateCellChat(ovary_young)
  ovary_old <- updateCellChat(ovary_old)
  object.list <- list(young = ovary_young, old = ovary_old)
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  
  
  #> Merge the following slots: 'data.signaling','images','net', 'netP','meta', 
  #> 'idents', 'var.features' , 'DB', and 'LR'.
  cellchat
  # An object of class CellChat created from a merged object with multiple datasets 
  # 649 signaling genes.
  # 2867 cells. 
  # CellChat analysis of single cell RNA-seq data!
  
  
  # Create a directory to save figures----
  data.dir <- './comparison'
  dir.create(data.dir)
  setwd(data.dir)
  
  # Part I: Predict general principles of cell-cell communication  有用 ------
  # 比较总体的相互作用数量及强度 - 柱状图
  gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  gg1 + gg2
  
  # 显示整体的细胞通讯数量及强度
  # Compare the number of interactions and interaction strength among different cell populations
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  # colored edges represent signaling in the second dataset compared to the first one
  gg1 <- netVisual_heatmap(cellchat)
  #> Do heatmap based on a merged object
  gg2 <- netVisual_heatmap(cellchat, measure = "weight")
  #> Do heatmap based on a merged object
  gg1 + gg2
  
  
  levels(cellchat@meta[["labels"]])
  # SC GC SMC Oocyte NKT EC Macrophage
  # ovary_young@net$count
  #              SC GC SMC Oocyte NKT  EC Macrophage
  # SC          51 11  33    111  49  76         33
  # GC          19 22  19     89  23  43         26
  # SMC         40 11  34     98  39  71         22
  # Oocyte     114 83 104    165  96 136        101
  # NKT         40 19  30     71  34  57         45
  # EC          77 34  56    117  58  89         56
  # Macrophage  43 32  30     96  57  75         43
  
  # 比较不同细胞类型之间的细胞通讯数量及强度
  weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_circle(object.list[[i]]@net$count, weight.scale = T, 
                     label.edge= T, # T则显示相应的Interations数量
                     edge.weight.max = weight.max[2], edge.width.max = 12, 
                     title.name = paste0("Number of interactions - ", names(object.list)[i]),
                     sources.use = c(4),
                     targets.use = c(1:7))
  }
  # 比较主要的sources 和 targets 
  num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
  gg <- list()
  for (i in 1:length(object.list)) {
    gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
  }
  patchwork::wrap_plots(plots = gg)
  ## Identify signaling changes associated with one cell group -------------------
  gg1 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "Oocyte", signaling.exclude = "MIF")
  gg2 <- netAnalysis_signalingChanges_scatter(cellchat, idents.use = "GC", signaling.exclude = c("MIF"))
  patchwork::wrap_plots(plots = list(gg1,gg2))
  
  # 比较ovary_young和ovary_old之间的AMH通路
  pathways.show <- c("AMH") # AMH,MK,PTN  ,因为这里是画两种细胞类型的比较图，所以选择的通路要在两者中均存在
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show, layout = "circle", 
                        edge.weight.max = weight.max[1], edge.width.max = 10, scale = TRUE,
                        signaling.name = paste(pathways.show, names(object.list)[i]))
  }
  
  netVisual_aggregate(ovary_young, signaling = "AMH", layout = "circle", scale = TRUE)
  
  
  # saveRDS(cellchat, file = "cellchat_comparisonAnalysis_ovary_young_vs_old.rds")
  cellchat <- readRDS("cellchat_comparisonAnalysis_ovary_young_vs_old.rds")