
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


LRpairs_list <- 'IL13_IL13RA1'
LRpairs <- as.character(unlist(strsplit(LRpairs_list,split=",")))
vertex.receiver=c(20,19,12,13,8,5,21,6,15,16)
vertex.receiver2=c(11,9,10,14,3,4,17,1,18,7,2)
plot_hierarchy(cellchat=cellchat,LRpairs=LRpairs,pathway_name='IL4',group.new=NULL,vertex.receiver=vertex.receiver,vertex.receiver2=vertex.receiver2)


###Figs12d
library(dplyr)
library(ggplot2)
library(readxl)

bar.df <- read_excel("metadata_all_blister.xlsx")

bar.df <- mutate(bar.df,name=factor(bar.df$Sample, levels=c("BP10","BP9","BP8","BP7")))
text.df <- as.data.frame(table(bar.df$name))
color_cluster=c("#FDBF6F","#1F78B4","#A6CEE3","#B2DF8A","#6A3D9A","#33A02C","#FFFF99","#B15928","#E31A1C","#CAB2D6","#FF7E00","#1B9E77","#FB9A99")
names(color_cluster)=c('Eosinophil','Effector CD8 T','Memory CD4 T','Treg','DC','Macrophage','Unknown','Proliferating T','Th2/Tc2','NK','pDC','Keratinocyte','Neutrophil')

ggplot(bar.df,aes(x=name))+
  geom_bar(aes(fill=CellType),position = "fill",width = .7)+
  scale_x_discrete("")+
  scale_y_continuous("Relative proportion of cells",expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("Cell types",values = color_cluster)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+
  coord_flip()


library(Seurat)
library(RColorBrewer)
cols=c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))
col=unique(cols)
col=col[-c(17,23)]
load('project.rdata')
pdf(file="UMAP_SC_BySample.pdf",width=20)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = TRUE,group.by = "Cluster",cols=col,split.by='Sample')+coord_fixed()
print(p)
dev.off()









library(CellChat)
load('cellchat.RData')
scource('functon.r')

LRpairs_list <- 'IL13_IL13RA1'
LRpairs <- as.character(unlist(strsplit(LRpairs_list,split=",")))
pathway_name <- 'IL4'
for(i in 1:length(LRpairs)){
    LRpairs_s <- data.frame(interaction_name=LRpairs[i])
    vertex.receiver=c(6,7,12,1,2,3,9)
    vertex.receiver2=c(10,4,8,5,11,13)
    pdf(paste0(LRpairs[i],"_netVisual_hierarchy.pdf"),width=15,height=15)
    netVisual_individual(cellchat, signaling=pathway_name,pairLR.use=LRpairs_s,vertex.receiver=vertex.receiver,vertex.receiver2=vertex.receiver2,layout = "hierarchy")
    dev.off()
}
