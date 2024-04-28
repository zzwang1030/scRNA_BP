library(Seurat)
library(ggplot2)
load('project.rdata')

barcode_info <- data.frame(Bracode=rownames(project@meta.data),project@meta.data[,c('Sample','Group','Cluster','CellType')])
write.table(barcode_info,'metadata_skin_immune_cells.xls',row.names=FALSE,sep='\t',quote=FALSE)
# umap
project$CellType <- factor(project$CellType,levels=c("Treg","Memory CD4 T","Effector CD8 T","Th2","Proliferating T","NK and gdT cells ","Langerhans cells","GPNMBhigh Macro","CD163+ Macro","FCN1+ Macro","CXCL1high DC","cDC2","LAMP3+ DC","pDC and Mast"))
col1<-c('#F8766D','#E38900','#C49A00','#99A800','#53B400','#00BC56','#00C094','#00BFC4','#00B6EB','#06A4FF','#A58AFF','#DF70F8','#FB61D7','#FF66A8')
col2 <- c('#FF8000','#008000')
pdf(file="UMAP_SC_ByCellType.pdf",width=10)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "CellType",cols=col1)
print(p)
dev.off()

pdf(file="UMAP_SC_ByGroup.pdf",width=7)
p <- DimPlot(object = project, reduction = "umap",label.size = 4, pt.size=0.5,label = FALSE,group.by = "Group",cols=col2)+coord_fixed()
print(p)
dev.off()
