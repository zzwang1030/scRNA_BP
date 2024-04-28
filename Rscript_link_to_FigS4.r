library(ggplot2)
data <- read.table("GO_enrichment_significant_topGO.xls", sep="\t", header=TRUE)
fdenrich <- data[, "Significant"] / data[, "Annotated"]
pvs <- data$P_value
data <- data[order(pvs),]
if(dim(data)[1] > 20){ data <- data[1:20,] }
pathname <- as.vector(data[, "GO_term"])
fdenrich <- data[, "Significant"] / data[, "Annotated"]
data <- data[order(fdenrich, decreasing=TRUE),]
pvs <- data[, "P_value"]
pvs <- -log10(pvs)
xmin <- min(fdenrich, na.rm=TRUE); xmin <- round(xmin, 3); xmax <- max(fdenrich, na.rm=TRUE); xmax <- round(xmax, 3); ymin <- min(pvs); ymin <- round(ymin, 0); ymax <- max(pvs); ymax <- round(ymax, 0)
ncolors <- rainbow(length(pathname))
ncounts <- data[, "Significant"]
aa <- data.frame(fdenrich=fdenrich, pathname=pathname, pvs=pvs, ncounts=ncounts)
aa <- aa[order(aa$fdenrich, decreasing = TRUE),] 
aa$pathname <- factor(aa$pathname, levels=rev(as.vector(aa$pathname)))
GeneCount <- aa$ncounts
Minuslog10Pvalue <- aa$pvs
p_2D <- ggplot(data=aa, aes(x=fdenrich, y=pathname))+geom_point(aes(colour=Minuslog10Pvalue, size = GeneCount), stat="identity", position="identity") + scale_colour_gradient(name='-log10(Pvalue)',low = "#92E78D", high="#4A4FAC")+scale_size(name='Gene Number')+theme_bw()+xlab("Rich Factor")+ylab("GO Term")+
theme(text=element_text(family="serif"),panel.grid=element_blank(), plot.title = element_text(hjust = 0.5, face="bold", colour="black", size=25),legend.text=element_text(size=20), axis.title.x=element_text(size=30),axis.title.y=element_text(size=30),axis.text.y=element_text(size=20),axis.text.x=element_text(size=20,colour="black", angle=90, hjust=1),legend.title = element_text(size = 30))+ggtitle("GO Enrichment Scatter Plot")
pdf(file="GO_enrichment_topGO_2D.pdf",height=15,width=18)
print(p_2D)
dev.off()
