
library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)


# Fig5a
load('cellchat.RData')
LRpairs <- 'CCL17_CCR4'
pathway_name <- cellchat@LR$LRsig[cellchat@LR$LRsig$interaction_name%in%LRpairs,]$pathway_name
LRpairs_s <- data.frame(interaction_name=LRpairs)
vertex.receiver=c(7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
pdf(paste0('All_',LRpairs,"_netVisual_hierarchy.pdf"),width=15,height=15)
netVisual_individual(cellchat, signaling=pathway_name,pairLR.use=LRpairs_s,layout = "hierarchy",vertex.receiver=vertex.receiver)
dev.off()