library(Seurat)
library(ggplot2)
load('project.rdata')

barcode_info <- data.frame(Bracode=rownames(project@meta.data),project@meta.data[,c('Sample','Group','CellType')])
write.table(barcode_info,'metadata.xls',row.names=FALSE,sep='\t',quote=FALSE)
# umap
project <- subset(project,idents=c(0,1,2,3,4,5,6,7,8,9))
project$CellType <- factor(project$CellType,levels=c("0_CCL19+ FB","1_APCDD1+ FB","2_CFH+ FB","3_HSPA1A+ FB","4_POSTN+ FB","5_ASPN+ FB","6_FBLN1+ FB","7_IGFBP2+ FB","8_INHBA+ FB","9_TAGLN+ FB"))
col1<-c('#4E79A7','#A0CBE8','#F28E2B','#FFBE7D','#59A14F','#8CD17D','#B6992D','#F1CE63','#499894','#95C4BF')
col2 <- c('#FF8000','#008000')
pdf(file="UMAP_SC_ByCellType.pdf")
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = TRUE,group.by = "CellType",cols=col1)+coord_fixed()
print(p)
dev.off()
pdf(file="UMAP_SC_BySample.pdf",width=40)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "CellType",split.by='Sample',cols=col1)+coord_fixed()
print(p)
dev.off()
pdf(file="UMAP_SC_ByGroup.pdf",width=7)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "Group",cols=col2)+coord_fixed()
print(p)
dev.off()


library(dplyr)
library(ggplot2)
library(readxl)

bar.df <- read_excel("metadata_skin_FB.xlsx")

bar.df <- mutate(bar.df,name=factor(bar.df$Sample, levels=c("HC8","HC7","HC6","HC5","HC4","HC3","HC2","HC1","BP5","BP4","BP3","BP2","BP1")))
text.df <- as.data.frame(table(bar.df$name))
color_cluster=c('#4E79A7','#A0CBE8','#F28E2B','#FFBE7D','#59A14F','#8CD17D','#B6992D','#F1CE63','#499894','#95C4BF')
names(color_cluster)=c("0_CCL19+ FB","1_APCDD1+ FB","2_CFH+ FB","3_HSPA1A+ FB","4_POSTN+ FB","5_ASPN+ FB","6_FBLN1+ FB","7_IGFBP2+ FB","8_INHBA+ FB","9_TAGLN+ FB")

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

bar.df <- read_excel("metadata_skin_KC.xlsx")

bar.df <- mutate(bar.df,name=factor(bar.df$Sample, levels=c("HC8","HC7","HC6","HC5","HC4","HC3","HC2","HC1","BP5","BP4","BP3","BP2","BP1")))
text.df <- as.data.frame(table(bar.df$name))
color_cluster=c('#F8766D','#DE8C00','#B79F00','#7CAE00','#00BA38','#00C08B','#00BFC4','#00B4F0',"#619CFF","#C77CFF","#F564E3","#FF64B0")
names(color_cluster)=c("0_Suprabasal","1_Basal","2_Suprabasal","3_Proliferating","4_Channels","5_Outer root sheath","6_Basal","7_Suprabasal","8_Inner root sheath","9_Stratum corneum","10_Stratum corneum","11_Basal")

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
library(ggplot2)
load('project.rdata')

barcode_info <- data.frame(Bracode=rownames(project@meta.data),project@meta.data[,c('Sample','Group','Cluster','CellType')])
write.table(barcode_info,'metadata.xls',row.names=FALSE,sep='\t',quote=FALSE)
# umap
col2 <- c('#FF8000','#008000')
pdf(file="UMAP_SC_ByCellType.pdf",width=8)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = TRUE,group.by = "CellType")
print(p)
dev.off()
pdf(file="UMAP_SC_ByGroup.pdf",width=7)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "Group",cols=col2)+coord_fixed()
print(p)
dev.off()