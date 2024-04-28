
library(Seurat)
library(ggplot2)
col1 <- c('#FF8000','#008000')
# a_dotplot
load('project.rdata')
marker <- read.table('All_PBMC_marker.txt',header=FALSE)
markers <- marker$V1
project@active.ident <- factor(project@active.ident,levels=c(9,1,3,8,7,4,2,12,5,13,10,11,6,14,15))
pdf('dotplot.pdf',height=5,width=17)
p <- DotPlot(project, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()
print(p)
dev.off()


# b_dotplot
load('project.rdata')
marker <- read.table('PBMC_Mermory_CD4_T_marker.txt',header=FALSE)
markers <- marker$V1
pdf('dotplot.pdf',height=2,width=7)
p <- DotPlot(project, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()
print(p)
dev.off()

# c_dotplot
load('project.rdata')
marker <- read.table('PBMC_Myeiold_marker.txt',header=FALSE)
markers <- marker$V1
pdf('dotplot.pdf',height=3,width=10)
p <- DotPlot(project, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()
print(p)
dev.off()

library(dplyr)
library(ggplot2)
library(readxl)

bar.df <- read_excel("metadata_all_pbmc.xlsx")

bar.df <- mutate(bar.df,name=factor(bar.df$Sample, levels=c("HC16","HC15","HC14","HC13","HC12","HC11","HC10","HC9","BP8","BP7","BP6","BP5","BP4","BP3","BP2","BP1")))
text.df <- as.data.frame(table(bar.df$name))
color_cluster=c("#FF7A9A","#05C49E","#00C6CE","#B899FF","#00C3F1","#00C5A7","#00BEFF","#75A9FF","#E29D3F","#00B5FF","#FF8075","#FB8E56","#C1AB36","#6CBD5B","#9BB642")
names(color_cluster)=c("NK-1","Memory CD4 T","Naive CD4 T","NK-2","B cells","Effector CD8 T","Monocytes","MAIT cells","Proliferating T","Naive CD8 T","Megakaryocytes","DC","pDC","unknown","Neutrophils")

ggplot(bar.df,aes(x=name))+
  geom_bar(aes(fill=CellType),position = "fill",width = .7)+
  scale_x_discrete("")+
  scale_y_continuous("Proportion of each cluster among all PBMC cells",expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("Cell types",values = color_cluster)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+
  coord_flip()


library(Seurat)
library(ggplot2)
load('project.rdata')

barcode_info <- data.frame(Bracode=rownames(project@meta.data),project@meta.data[,c('Sample','Group','Cluster','CellType')])
write.table(barcode_info,'metadata.xls',row.names=FALSE,sep='\t',quote=FALSE)
# umap
col2 <- c('#FF8000','#008000')
pdf(file="UMAP_SC_BySample.pdf",width=8)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = TRUE,group.by = "Sample")
print(p)
dev.off()
