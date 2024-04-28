
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)


source('netVisual_individual.R')
load('1.1_cellchat/Control/cellchat.RData')
cellchat_1 <- cellchat 
load('1.1_cellchat/BP/cellchat.RData')
cellchat_2 <- cellchat
LRpairs <- 'CXCL12_CXCR4'
pathway_name <- cellchat_1@LR$LRsig[cellchat_1@LR$LRsig$interaction_name%in%LRpairs,]$pathway_name
LRpairs_s <- data.frame(interaction_name=LRpairs)
vertex.receiver=c(22,23)
pdf(paste0('Control_',LRpairs,"_netVisual_hierarchy.pdf"),width=15,height=15)
netVisual_individual(cellchat_1, signaling=pathway_name,pairLR.use=LRpairs_s,layout = "hierarchy",vertex.receiver=vertex.receiver)
dev.off()
pathway_name <- cellchat_2@LR$LRsig[cellchat_2@LR$LRsig$interaction_name%in%LRpairs,]$pathway_name
LRpairs_s <- data.frame(interaction_name=LRpairs)
vertex.receiver=c(22,23)
pdf(paste0('BP_',LRpairs,"_netVisual_hierarchy.pdf"),width=15,height=15)
netVisual_individual(cellchat_2, signaling=pathway_name,pairLR.use=LRpairs_s,layout = "hierarchy",vertex.receiver=vertex.receiver)
dev.off()
Sample1Name <- 'Control'
Sample2Name <- 'BP'
levels(cellchat_1@idents)
levels(cellchat_2@idents)
types <- intersect(levels(cellchat_1@idents),levels(cellchat_2@idents))
object.list <- list(S1=cellchat_1,S2=cellchat_2)
names(object.list) <- c(Sample1Name,Sample2Name)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
LRpairs_filter <- data.frame(interaction_name='CXCL12_CXCR4')
pathway_name <- object.list[[1]]@LR$LRsig[object.list[[1]]@LR$LRsig$interaction_name%in%'CXCL12_CXCR4',]$pathway_name
for (i in 1:length(object.list)) {
  pdf(paste0(names(object.list)[i],"_CXCL12_CXCR4_netVisual_circle.pdf"),width=7,height=7)
  netVisual_individual(object.list[[i]],signaling=pathway_name, pairLR.use=LRpairs_filter,layout = "circle",signaling.name = names(object.list)[i])
  dev.off()
} 


###Fig4i
library(Seurat)
library(ggplot2)
load('/2_QC_DEG_plot_SubFibroblasts/project.rdata')
project_PLA2G2A <- subset(project, subset = PLA2G2A > 0)
project_CXCL12 <- subset(project, subset = CXCL12 > 0)
project_PLA2G2A_CXCL12 <- subset(project, subset = PLA2G2A > 0 & CXCL12 > 0)
project_DN <- subset(project, subset = PLA2G2A == 0 & CXCL12 == 0)
PLA2G2A <- paste0('PLA2G2A+ ',dim(project_PLA2G2A)[2]-dim(project_PLA2G2A_CXCL12)[2],'cells')
CXCL12 <- paste0('CXCL12+ ',dim(project_CXCL12)[2]-dim(project_PLA2G2A_CXCL12)[2],'cells')
DN <- paste0('DN ',dim(project_DN)[2],'cells')
PLA2G2A_CXCL12 <- paste0('PLA2G2A+CXCL12+ ',dim(project_PLA2G2A_CXCL12)[2],'cells')
project$CellType <- DN
project@meta.data[colnames(project_PLA2G2A_CXCL12),]$CellType <- PLA2G2A_CXCL12
project@meta.data[setdiff(colnames(project_PLA2G2A),colnames(project_PLA2G2A_CXCL12)),]$CellType <- PLA2G2A
project@meta.data[setdiff(colnames(project_CXCL12),colnames(project_PLA2G2A_CXCL12)),]$CellType <- CXCL12
project$CellType <- factor(project$CellType,levels=c(DN,PLA2G2A,CXCL12,PLA2G2A_CXCL12))
project.list <- SplitObject(object = project, split.by = "Group")
cols<-c('#BEBEBE','#00FF00','#FF0000','#0000FF')
pdf(file="UMAP.pdf",width=23,height=5)
p1 <- DimPlot(object = project,reduction = "umap",label.size = 4,pt.size=0.5,label = FALSE,group.by = "CellType",cols =cols)+ggtitle('all FB cells')
p2 <- DimPlot(object = project.list[['BP']],reduction = "umap",label.size = 4,pt.size=0.5,label = FALSE,cols =cols,group.by='CellType')+ggtitle('BP')
p3 <- DimPlot(object = project.list[['Control']],reduction = "umap",label.size = 4,pt.size=0.5,label = FALSE,cols =cols,group.by='CellType')+ggtitle('Control')
p<-plot_grid(p1,p2,p3,nrow=1)
print(p)
dev.off()
