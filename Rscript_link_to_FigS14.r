library(CellChat)
library(patchwork)
library(ComplexHeatmap)
library(dplyr)

##Figs14a
load('skin/All/cellchat.RData')
source('/netVisual_individual.R')
LRpairs <- c('CCL19_CCR7')
for(i in 1:length(LRpairs)){
  print(i)
  pathway_name <- cellchat@LR$LRsig[cellchat@LR$LRsig$interaction_name%in%LRpairs[i],]$pathway_name
  LRpairs_s <- data.frame(interaction_name=LRpairs[i])
  if(i==1){
    vertex.receiver=c(22,23)
  }else if(i==2){
    vertex.receiver=c(7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
  }
  pdf(paste0(LRpairs[i],"_netVisual_hierarchy.pdf"),width=15,height=15)
  netVisual_individual(cellchat, signaling=pathway_name,pairLR.use=LRpairs_s,layout = "hierarchy",vertex.receiver=vertex.receiver)
  dev.off()
}

##Figs14b
library(CellChat)
load('cellchat.RData')
scource('functon.r')
LRpairs_list <- 'CCL17_CCR4,CCL19_CCR7'
LRpairs <- as.character(unlist(strsplit(LRpairs_list,split=",")))
pathway_name <- 'CCL'
for(i in 1:length(LRpairs)){
    LRpairs_s <- data.frame(interaction_name=LRpairs[i])
    vertex.receiver=c(10,4,8,5,11,13,7)
    vertex.receiver2=c(6,12,1,2,3,9)
    pdf(paste0(LRpairs[i],"_netVisual_hierarchy.pdf"),width=15,height=15)
    netVisual_individual(cellchat, signaling=pathway_name,pairLR.use=LRpairs_s,vertex.receiver=vertex.receiver,vertex.receiver2=vertex.receiver2,layout = "hierarchy")
    dev.off()
}

