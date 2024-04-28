##Quality Control rCASC mitoRiboUmi ##
library("argparser")
library("refGenome")
library('Seurat')
library('data.table')
gtf.name<-'genes.gtf'
biotype <- 'protein_coding'
umiXgene<-'3'
mainMatrix <- read.table('gene_cell_exprs_table_symbol.xls',header=TRUE,sep='\t',row.names=1,check.names = FALSE)
beg <- ensemblGenome()
basedir(beg) <- getwd()
read.gtf(beg, gtf.name)
annotation <- extractPaGenes(beg)
rps <- grep("^RPS",toupper(rownames(mainMatrix)))
rpl <- grep("^RPL",toupper(rownames(mainMatrix)))
ribosomal <- c(rps, rpl)
ribosomal=rownames(mainMatrix)[ribosomal]
mitocondrial=intersect(annotation$gene_name[which(annotation$seqid=="MT")],rownames(mainMatrix))
mainMatrix=mainMatrix[union(union(intersect(annotation$gene_name[which(annotation$gene_biotype%in%biotype)],rownames(mainMatrix)),ribosomal),mitocondrial),]
x=colSums(mainMatrix[ribosomal,])/colSums(mainMatrix)*100
y=colSums(mainMatrix[mitocondrial,])/colSums(mainMatrix)*100
write.table(data.frame(Barcode=names(x), percent_ribosomal_count=x,percent_mitochondrial_count=y) ,'Barcode_ribomitoumi.xls',row.names=FALSE,col.names=TRUE,sep='\t',quote=FALSE)

#####Run normalition, dimensionality reduction, FindClusters and so on
library(Seurat)
library(data.table)
library(harmony)
mit_file <- 'Barcode_ribomitoumi.xls' # rCASC itoRiboUmi result
double_file <- 'predicted_label.xls' # Scrublet result
gene_matrix <- 'gene_cell_exprs_table_symbol.xls' # Matrix
sample_list <- 'BP1_skin,BP3_skin,BP2_skin,BP4_skin,BP5_skin,HC8_skin,HC7_skin,HC6_skin,HC4_skin,HC5_skin,HC3_skin,HC2_skin,HC1_skin'
#data import
project_mx <- fread(gene_matrix,sep="\t", header=TRUE,data.table=F)
colnames(project_mx) <- gsub("\\.", "-", colnames(project_mx))
rownames(project_mx) <- project_mx[,1]
project_mx[,1] <- NULL
project <- CreateSeuratObject(counts = project_mx, min.cells = 1, min.features = 0, project = "project")
project[["percent.redcell"]] <- PercentageFeatureSet(project, features=c('HBA1','HBA2','HBB'))
double1 <- read.table(double_file,header=TRUE,sep="\t",row.names=1)
barcode_info <- cbind(barcode_info,double1[colnames(project),])
index <- which(double1$predicted_doublets %in% 0)
valid_cells <- colnames(project[,which(colnames(project) %in% rownames(double1[index,]))])
mit_data <- read.table(mit_file,header=TRUE)
rownames(mit_data) <- mit_data$Barcode
mit_data <- mit_data[colnames(project),]
project[['percent_mitochondrial_count']] <- mit_data$percent_mitochondrial_count
samples <- as.character(unlist(strsplit(sample_list,split=",")))
project[["Sample"]]<-rep(0,dim(project)[2])
for(i in 1:length(samples)){
    yy<-paste("-",i,"$",sep="")
    project[["Sample"]][grep(yy,colnames(project)),]<-samples[i]
}
#data filter
project <- subset(project,cells=valid_cells)
project <- subset(project, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent_mitochondrial_count < 10 & percent.redcell < 10)
#harmony,dimensionality reduction, clustering, DEGs
DefaultAssay(project) <- 'RNA'
project <- NormalizeData(object = project, scale.factor = 10000)
project <- FindVariableFeatures(object = project,selection.method = "vst", nfeatures = 2000, mean.cutoff = c(0.1, Inf),dispersion.cutoff=c(1, Inf))
project <- ScaleData(object = project, verbose = TRUE,features=rownames(project[['RNA']]@data),assay='RNA')
project <- RunPCA(object = project, npcs = 30, verbose = TRUE,assay='RNA',features=VariableFeatures(project,assay='RNA'))
project <- RunHarmony(project, "Sample")
project <- FindNeighbors(object = project, reduction = "harmony", dims = 1:30,assay='RNA')
project <- RunTSNE(object = project, reduction = "harmony",dims=1:30,assay='RNA',check_duplicates=FALSE)
project <- RunUMAP(object = project, reduction = "harmony",dims = 1:30,assay='RNA')
project <- FindClusters(project,resolution = 0.2)
diffs <- FindAllMarkers(object = project, logfc.threshold = 0.1, only.pos = args$onlypos, test.use="bimod", min.pct = 0.01,assay="RNA")
save(list=c("project","diffs"),file='project.rdata')
              
