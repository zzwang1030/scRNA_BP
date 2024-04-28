
load('project.rdata')
#umap
temp <- c('#FF8000','#008000')
dev.off()
pdf(file="UMAP_SC_ByGroup.pdf")
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "Group",cols=temp)+coord_fixed()
print(p)
dev.off()
pdf(file="UMAP_SC_ByCluster.pdf")
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = TRUE,group.by = "Cluster")+coord_fixed()
print(p)
dev.off()
pdf(file="UMAP_SC_BySample.pdf",height=4,width=30)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,split.by='Sample',group.by = "Cluster")+coord_fixed()
print(p)
dev.off()
#dotplot
marker_list <- 'CD93, ACKR1, AQP1, TYRP1, PMEL, DCT, MLANA, CRABP1, CCL21, LYVE1, TFF3, TFPI, MMRN1, TAGLN, ACTA2, TPM2, MYL9, DCN, COL1A1, COL1A2, COL3A1, CFD, APOD, DCD, KRT19, AQP5, KRT7, SNORC, SCGB1B2P, SCGB1D2, CD68, CD1C, CLEC10A, LYZ, TYROBP, CD207, CXCL8, CD3D, CD3G, CD8A, PTPRC, GNLY, NKG7, FLG, LOR, KRT1, KRT10, S100A7, KRT6A, GJB6, KRT5, KRT14'
markers <- as.character(unlist(strsplit(marker_list,split=", ")))
height <- length(levels(project))*0.2+2
width <- length(markers)*0.25+1
pdf(file='bubble_plot.pdf',height=height,width=width)
p <- DotPlot(project, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()+
		theme(axis.text=element_text(size=8),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=8),legend.title=element_text(size=8),legend.key.size=unit(0.13, "inches"))+
		labs(y='CellType')
print(p)
dev.off()


library(dplyr)
library(ggplot2)
library(readxl)
bar.df <- read_excel("metadata_all_skin.xlsx")

bar.df <- mutate(bar.df,name=factor(bar.df$Sample, levels=c("HC8","HC7","HC6","HC5","HC4","HC3","HC2","HC1","BP5","BP4","BP3","BP2","BP1")))
text.df <- as.data.frame(table(bar.df$name))
color_cluster=c("#f8766d","#7CAE00","#7CAE00","#00BF7D","#A3A500","#00B0F6","#00BAE0","#C77CFF","#FF62BC")
names(color_cluster)=c("Keratinocytes","Fibroblasts","T/NK cells","Endothelial cells","Smooth muscle cells","DC/Mac cells","Melanocytes","Lymphatics cells","Sweet gland cells")

ggplot(bar.df,aes(x=name))+
  geom_bar(aes(fill=Celltype_annotation),position = "fill",width = .7)+
  scale_x_discrete("")+
  scale_y_continuous("Relative proportion of cells",expand = c(0,0),labels = scales::label_percent(),position = "right")+
  scale_fill_manual("Cell types",values = color_cluster)+
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent",colour = NA),
    axis.line.x = element_line(colour = "black")
  )+
  coord_flip() 