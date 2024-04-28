library(Seurat)
library(ggplot2)
load('project.rdata')

barcode_info <- data.frame(Bracode=rownames(project@meta.data),project@meta.data[,c('Sample','Group','Cluster','CellType')])
write.table(barcode_info,'metadata.xls',row.names=FALSE,sep='\t',quote=FALSE)
# umap
project$CellType <- factor(project$CellType,levels=c("Treg","Memory CD4 T","Effector CD8 T","Th2","Proliferating T","NK and gdT cells ","Langerhans cells","GPNMBhigh Macro","CD163+ Macro","FCN1+ Macro","CXCL1high DC","cDC2","LAMP3+ DC","pDC and Mast"))
col1<-c('#F8766D','#E38900','#C49A00','#99A800','#53B400','#00BC56','#00C094','#00BFC4','#00B6EB','#06A4FF','#A58AFF','#DF70F8','#FB61D7','#FF66A8')
col2 <- c('#FF8000','#008000')

pdf(file="UMAP_SC_BySample.pdf",width=40)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "CellType",split.by='Sample',cols=col1)+coord_fixed()
print(p)

# dotplot
celltype_file <- read.csv('SKIN_SubT_Myeloid_resolution2.csv',header=TRUE)
valid_cells <- intersect(celltype_file$Barcode,colnames(project))
celltype_file <- celltype_file[celltype_file$Barcode%in%valid_cells,]
rownames(celltype_file) <- celltype_file$Barcode
marker_list <- 'CD3D, CD3E, IL7R, CXCR4, CD8A, CD8B, IFNG, CCL5, GZMA, NKG7, XCL2, XCL1, AREG, TRDC, FOXP3, IL2RA, CTLA4, TNFRSF4, GATA3, IL17RB, IL13, IL5, MKI67, STMN1'
markers <- as.character(unlist(strsplit(marker_list,split=", ")))
project_sob <- subset(project,cells=celltype_file[celltype_file$Bubble%in%'T',]$Barcode)
celltype_file_sob <- celltype_file[colnames(project_sob),]
project_sob$CellType <- factor(celltype_file_sob$Celltype,levels=c('Memory CD4 T','Effector CD8 T','NK and gdT','Treg','Th2','Proliferating T'))
project_sob@active.ident <- project_sob$CellType
height <- length(levels(project_sob))*0.2+2
width <- length(markers)*0.25+1
pdf(file='FigS2a_dotplot.pdf',height=height,width=width)
p <- DotPlot(project_sob, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()+
		theme(axis.text=element_text(size=8),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=8),legend.title=element_text(size=8),legend.key.size=unit(0.13, "inches"))+
		labs(y='CellType')
print(p)
dev.off()

marker_list <- 'JCHAIN, CPA3, TCF4, GZMB, TPSB2, TPSAB1, CD1A, CD207, HLA-DQB2, FCER1A, CLEC10A, CD1C, LAMP3, CCL17, CCL22, MMP12, CXCL1, CXCL8, CXCL3, CD68, CD14, FCGR3A, CD163, CCL18, CCL13, VCAN, FCN1, S100A12, S100A9, S100A8, APOE, CD9, GPNMB, S100A2, CCL27'
markers <- as.character(unlist(strsplit(marker_list,split=", ")))
project_sob <- subset(project,cells=celltype_file[celltype_file$Bubble%in%'Monocyte',]$Barcode)
celltype_file_sob <- celltype_file[colnames(project_sob),]
project_sob$CellType <- factor(celltype_file_sob$Celltype,levels=c('pDC and Mast','Langerhans cells','cDC2','LAMP3+ DC','CXCL1high DC','CD163+ Macro','FCN1+ Macro','GPNMBhigh Macro'))
project_sob@active.ident <- project_sob$CellType
height <- length(levels(project_sob))*0.2+2
width <- length(markers)*0.25+1
pdf(file='FigS2b_dotplot.pdf',height=height,width=width)
p <- DotPlot(project_sob, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()+
		theme(axis.text=element_text(size=8),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=8),legend.title=element_text(size=8),legend.key.size=unit(0.13, "inches"))+
		labs(y='CellType')
print(p)
dev.off()




library(dplyr)
library(ggplot2)
library(readxl)

bar.df <- read_excel("metadata_skin_immune_cells_T.xlsx")

bar.df <- mutate(bar.df,name=factor(bar.df$Sample, levels=c("HC8","HC7","HC6","HC5","HC4","HC3","HC2","HC1","BP5","BP4","BP3","BP2","BP1")))
text.df <- as.data.frame(table(bar.df$name))
color_cluster=c("#F8766D","#E38900","#C49A00","#99A800","#53B400","#00BC56")
names(color_cluster)=c("Treg","Memory CD4 T","Effector CD8 T","Th2","Proliferating T","NK and gdT cells")

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


library(dplyr)
library(ggplot2)
library(readxl)

bar.df <- read_excel("metadata_skin_immune_cells_M.xlsx")

bar.df <- mutate(bar.df,name=factor(bar.df$Sample, levels=c("HC8","HC7","HC6","HC5","HC4","HC3","HC2","HC1","BP5","BP4","BP3","BP2","BP1")))
text.df <- as.data.frame(table(bar.df$name))
color_cluster=c('#00C094','#00BFC4','#00B6EB','#06A4FF','#A58AFF','#DF70F8','#FB61D7','#FF66A8')
names(color_cluster)=c("Langerhans cells","GPNMBhigh Macro","CD163+ Macro","FCN1+ Macro","CXCL1high DC","cDC2","LAMP3+ DC","pDC and Mast")

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