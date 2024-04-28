##FigS13a
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

LRpairs_list <- 'CCL17_CCR4'
LRpairs <- as.character(unlist(strsplit(LRpairs_list,split=",")))
vertex.receiver=c(11,9,10,14,3,4,17,1,18,7,2)
vertex.receiver2=c(20,19,12,13,8,5,21,6,15,16)
plot_hierarchy(cellchat=cellchat,LRpairs=LRpairs,pathway_name='CCL',group.new=NULL,vertex.receiver=vertex.receiver,vertex.receiver2=vertex.receiver2)



##FigS13b
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


