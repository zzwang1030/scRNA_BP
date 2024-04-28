library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)

data_c <- read.table('/Control/CellChat/celltype_and_cell_interaction_result.xls',header=TRUE,sep='\t',quote='')
interationPair_c <- unique(data_c[data_c$annotation%in%'Secreted Signaling',]$interationPair)
data_b <- read.table('/BP/CellChat/celltype_and_cell_interaction_result.xls',header=TRUE,sep='\t',quote='')
interationPair_b <- unique(data_b[data_b$annotation%in%'Secreted Signaling',]$interationPair)
diff_b <- setdiff(interationPair_b,interationPair_c) 
inter <- intersect(interationPair_b,interationPair_c) 
all <-data.frame(interaction_name=interationPair_b)
load('/merge/cellchat_merge.RData')
pdf('Fibroblast0.pdf',width=12,height=35)
netVisual_bubble(cellchat, sources.use = 22,  comparison = c(2,1), angle.x = 45,pairLR.use=all)
dev.off()
pdf('Fibroblast1.pdf',width=12,height=35)
netVisual_bubble(cellchat, sources.use = 23,  comparison = c(2,1), angle.x = 45,pairLR.use=all)
dev.off()



