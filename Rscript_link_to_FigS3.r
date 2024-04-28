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

pdf(file="UMAP_SC_BySample.pdf",width=40)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "CellType",split.by='Sample',cols=col1)+coord_fixed()
print(p)
dev.off()

# dotplot # 
marker_list <- 'CCL19, CD74, C3, APOE, APOD, CFH, CFD, HSPA1A, HSPA1B, CXCL3, WIF1, APCDD1, COL18A1, COL13A1, COMP, COL1A1, POSTN, COL11A1, ASPN, OGN, CRABP1, IGFBP5, FBLN1, CCN5, SLPI, MFAP5, ANGPTL1, TIMP3, ENTPD1, IGFBP2, FGFBP2, SPON2, SERPINE2, RSPO3, INHBA, TAGLN, NR2F2, KLF5'
markers <- as.character(unlist(strsplit(marker_list,split=", ")))
project@active.ident <- factor(project@active.ident,levels=c("0_CCL19+ FB","2_CFH+ FB","3_HSPA1A+ FB","1_WIF1+ FB","4_POSTN+ FB","5_ASPN+ FB","6_FBLN1+ FB","7_IGFBP2+ FB","8_INHBA+ FB","9_TAGLN+ FB"))
height <- length(levels(project))*0.2+2
width <- length(markers)*0.25+1
pdf(file='FigS2c_dotplot.pdf',height=height,width=width)
p <- DotPlot(project, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()+
		theme(axis.text=element_text(size=8),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=8),legend.title=element_text(size=8),legend.key.size=unit(0.13, "inches"))+
		labs(y='CellType')
print(p)
dev.off()



library(Seurat)
library(ggplot2)
load('project.rdata')


barcode_info <- data.frame(Bracode=rownames(project@meta.data),project@meta.data[,c('Sample','Group','Cluster','CellType')])
write.table(barcode_info,'metadata.xls',row.names=FALSE,sep='\t',quote=FALSE)
# umap
col2 <- c('#FF8000','#008000')

pdf(file="UMAP_SC_BySample.pdf",width=40)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "CellType",split.by='Sample')+coord_fixed()
print(p)
dev.off()

# dotplot
marker_list<-'ASPRV1, LOR, FLG2, FLG, KRT2, KRT10, KRT1, S100A7, S100A9, SPRR1B, SERPINB4, S100A2, KRT15, KRT5, KRT14, COL17A1, SFRP1, CHI3L1, CCL2, MKI67, STMN1, UBE2C, TOP2A, KRT6B, KRT6C, APOC1, CD74, COL1A1, COL1A2, COL3A1, SAT1, GJB6, ATP1B1, ATP1A1'
markers <- as.character(unlist(strsplit(marker_list,split=", ")))
project@active.ident <- factor(project@active.ident,levels=c("9_Stratum corneum","10_Stratum corneum","0_Suprabasal","2_Suprabasal","7_Suprabasal","1_Basal","6_Basal","11_Basal","3_Proliferating","5_Outer root sheath","8_Inner root sheath","4_Channels"))
height <- length(levels(project))*0.2+2
width <- length(markers)*0.25+1
pdf(file='FigS2d_dotplot.pdf',height=height,width=width)
p <- DotPlot(project, features = as.character(markers),assay="RNA",cols = c("grey","red")) + RotatedAxis()+
		theme(axis.text=element_text(size=8),axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),legend.text=element_text(size=8),legend.title=element_text(size=8),legend.key.size=unit(0.13, "inches"))+
		labs(y='CellType')
print(p)
dev.off()














