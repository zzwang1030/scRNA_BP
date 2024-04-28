library(Seurat)

load('project.rdata')
sample_list <- 'BP1,BP2,BP3,BP4,BP5,BP6,BP7,BP8,HC9,HC10,HC11,HC12,HC13,HC14,HC15,HC16'
samples <- as.character(unlist(strsplit(sample_list,split=",")))
project[["Sample"]]<-rep(0,dim(project)[2])
for(i in 1:length(samples)){
    yy<-paste("-",i,"$",sep="")
    project[["Sample"]][grep(yy,colnames(project)),]<-samples[i]
}
project$Sample <-factor(project$Sample,levels=c("BP1","BP2","BP3","BP4","BP5","BP6","BP7","BP8","HC9","HC10","HC11","HC12","HC13","HC14","HC15","HC16"))
QC_index<-c("nFeature_RNA", "nCount_RNA", "percent_mitochondrial_count")
pdf(file="QC_index_by_Violinplot_Sample.pdf",width=5*length(QC_index)+2)
p<-VlnPlot(project, features = QC_index, ncol = length(QC_index),group.by = "Sample",pt.size=0)
print(p)
dev.off()



library(Seurat)
# 1_QC_vlnplot
load('project.rdata')
sample_list <- 'BP9,BP10,BP7,BP8'
samples <- as.character(unlist(strsplit(sample_list,split=",")))
project[["Sample"]]<-rep(0,dim(project)[2])
for(i in 1:length(samples)){
    yy<-paste("-",i,"$",sep="")
    project[["Sample"]][grep(yy,colnames(project)),]<-samples[i]
    
}
project$Sample <-factor(project$Sample,levels=c('BP7','BP8','BP9','BP10'))
QC_index<-c("nFeature_RNA", "nCount_RNA", "percent_mitochondrial_count")
pdf(file="QC_index_by_Violinplot_Sample.pdf",width=5*length(QC_index)+2)
p<-VlnPlot(project, features = QC_index, ncol = length(QC_index),group.by = "Sample",pt.size=0)
print(p)
dev.off()




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


