library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)

#Fig3a
load('/merge/cellchat_merge.RData')
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, return.data = T,color.use=c('#00BFC4','#F8766D'))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, return.data = T,color.use=c('#00BFC4','#F8766D'))
pdf("2.1_rankNet_Signaling_comparison.pdf",width = 14,height = 10);print(gg1[[2]] + gg2[[2]]);dev.off()


#Fig3b
legend_labels <- c('Th2','Treg','Memory CD4 T','Effector CD8 T','Proliferating T','NK and gdT','Langerhans_cells','GPNMBhigh Macro','CD163+ Macro','FCN1+ Macro','CXCL1high DC','cDC2','LAMP3+ DC','pDC and Mast','0_Suprabasal','1_Basal','3_Proliferating KC','5_Outer root sheath','6_Basal','7_Suprabasal','0_CCL19+ FB','1_APCDD1+ FB')
Sample1Name <- 'Control'
Sample2Name <- 'BP'
load('/Control/cellchat.RData')
cellchat_1 <- cellchat 
load('/BP/cellchat.RData')
cellchat_2 <- cellchat 
levels(cellchat_1@idents) 
levels(cellchat_2@idents)
object.list <- list(S1=cellchat_1,S2=cellchat_2)
names(object.list) <- c(Sample1Name,Sample2Name)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, return.data = T)
pathways.show <- as.character(unique(gg1$signaling.contribution$name))
path <- 'IL4'
pdf(paste0(path,"_netVisual_circle.pdf"),width=22,height=6)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = path, layout = "circle",  edge.width.max = 10, signaling.name = paste(path, names(object.list)[i]),thresh = 1)
    legend("right", legend=legend_labels, col=scPalette(length(legend_labels)),bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE,inset=c(0, 0.1))
}
dev.off()

# Fig3e
load('/BP/cellchat.RData')
LRpairs_filter <- cellchat@LR$LRsig[cellchat@LR$LRsig$pathway_name%in%'IL4',]$interaction_name
LRpairs_filter <- data.frame(interaction_name=LRpairs_filter)
pdf("BP_IL4_bubble_plot.pdf",width = 7,height = 3)
p <- netVisual_bubble(cellchat, sources.use = c(1),thresh = 1,targets.use = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),pairLR.use = LRpairs_filter, remove.isolate = F,line.on = T)
print(p)
dev.off()
load('/NO_527/3_redo_cellchat/1_test/1.1_cellchat/Control/cellchat.RData')
LRpairs_filter <- cellchat@LR$LRsig[cellchat@LR$LRsig$pathway_name%in%'IL4',]$interaction_name
LRpairs_filter <- data.frame(interaction_name=LRpairs_filter)
pdf("Control_IL4_bubble_plot.pdf",width = 7,height = 2.5)
p <- netVisual_bubble(cellchat, sources.use = c(1),thresh = 1,targets.use = c(2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23),pairLR.use = LRpairs_filter, remove.isolate = F,line.on = T)
print(p)
dev.off()

# Fig3f 
LRpairs <- 'IL13_IL13RA1'
load('/BP/cellchat.RData')
pathway_name <- cellchat@LR$LRsig[cellchat@LR$LRsig$interaction_name%in%LRpairs,]$pathway_name
LRpairs_s <- data.frame(interaction_name=LRpairs)
vertex.receiver=c(1,5,6,2,3,4)
pdf(paste0('BP_',LRpairs,"_netVisual_hierarchy.pdf"),width=15,height=15)
netVisual_individual(cellchat, signaling=pathway_name,pairLR.use=LRpairs_s,layout = "hierarchy",vertex.receiver=vertex.receiver)
dev.off()
load('/BP/cellchat.RData')
pathway_name <- cellchat@LR$LRsig[cellchat@LR$LRsig$interaction_name%in%LRpairs,]$pathway_name
LRpairs_s <- data.frame(interaction_name=LRpairs)
vertex.receiver=c(1,5,6,2,3,4)
pdf(paste0('Control_',LRpairs,"_netVisual_hierarchy.pdf"),width=15,height=15)
netVisual_individual(cellchat, signaling=pathway_name,pairLR.use=LRpairs_s,layout = "hierarchy",vertex.receiver=vertex.receiver)
dev.off()