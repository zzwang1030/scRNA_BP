##Fig6b
library(CellChat)
Sample1Name <- 'PBMC-ontrol'
Sample2Name <- 'PBMC-BP'
load('3.1_cellchat/PBMC-ontrol/cellchat.RData')
cellchat_1 <- cellchat 
load('3.1_cellchat/PBMC-BP/cellchat.RData')
cellchat_2 <- cellchat 
levels(cellchat_1@idents) 
levels(cellchat_2@idents)
object.list <- list(S1=cellchat_1,S2=cellchat_2)
names(object.list) <- c(Sample1Name,Sample2Name)
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = path)
pdf(paste0(path,"_netVisual_circle.pdf"),width=12,height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = path, layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(path, names(object.list)[i]))
}
dev.off()

##Fig6c
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)
scource('functon.r')
plot_hierarchy <- function(cellchat,LRpairs,pathway_name,group.new,vertex.receiver,vertex.receiver2){
    if(!is.null(group.new)){cellchat <- liftCellChat(cellchat, group.new)}
    LRpairs <- intersect(LRpairs,dimnames(cellchat@net$prob)[[3]])
    for(i in 1:length(LRpairs)){
        LRpairs_s <- data.frame(interaction_name=LRpairs[i])
        pdf(paste0(LRpairs[i],"_netVisual_hierarchy.pdf"),width=15,height=15)
        p <- netVisual_individual(cellchat, signaling=pathway_name,pairLR.use=LRpairs_s,vertex.receiver=vertex.receiver,vertex.receiver2=vertex.receiver2,layout = "hierarchy")
        print(p)
        dev.off()
    }
}
load('cellchat.RData')
path<-'IL4'
pdf(paste0(path,"_netVisual_bubble.pdf"),width=8,height=4)
netVisual_bubble(cellchat, sources.use = c(19,20), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,21), signaling = path, remove.isolate = FALSE)
dev.off()



##Fig6d
library(Seurat)
library(RColorBrewer)
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)
col=col[-c(17,23)]
load('project.rdata')
pdf(file="UMAP_SC_ByCluster.pdf",width=8)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = TRUE,group.by = "Cluster",cols=col)
print(p)
dev.off()


##Fig6e
library(Seurat)
library(ggplot2)
col1 <- c('#FF8000','#008000')
load('project.rdata')
marker <- read.table('blister_marker.txt',header=FALSE)
markers <- marker$V1
project@active.ident <- factor(project@active.ident,levels=c(1,2,9,3,12,6,7,5,11,4,10,8,13))
pdf('dotplot.pdf',height=5,width=13)
p <- DotPlot(project, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()
print(p)
dev.off()

##Fig6f and g
library(CellChat)
load('cellchat.RData')
scource('functon.r')
pathway_name<-'IL4'
pdf(paste0(pathway_name,"_netVisual_circle.pdf"),width=11,height=6)
netVisual_aggregate(cellchat, signaling = pathway_name, layout = "circle")
dev.off() 
pdf(paste0(pathway_name,"_netVisual_bubble.pdf"),width=8,height=4)
netVisual_bubble(cellchat, sources.use = c(6,7,12), targets.use = c(1,2,3,4,5,6,8,9,10,11,13), signaling = pathway_name, remove.isolate = FALSE)
dev.off()


