
library(reshape2)
library(ktplots)
trace('plot_cpdb', edit = T, where = asNamespace("ktplots"))
mypython <- '/software/Miniconda3/envs/cpdb/bin/cellphonedb'

### Th2
data <- read.table('/Th2.txt',header=TRUE,sep='\t',quote='')
data$CellType <- paste0(data$CellType_a,'|',data$CellType_b)
data$partner_a <- 'NULL'
data$partner_b <- 'NULL'
data$secreted <- 'NULL'
data$is_integrin <- 'NULL'
data$receptor_a <- 'NULL'
data$receptor_b <- 'NULL'
wide_mean <- dcast(data, id_cp_interaction+interacting_pair+partner_a+partner_b+gene_a+gene_b+secreted+receptor_a+receptor_b+annotation_strategy+is_integrin~ CellType, value.var = "mean")
wide_mean[is.na(wide_mean)] <- 0
wide_pvalue <- dcast(data, id_cp_interaction+interacting_pair+partner_a+partner_b+gene_a+gene_b+secreted+receptor_a+receptor_b+annotation_strategy+is_integrin~ CellType, value.var = "pValue")
wide_pvalue[is.na(wide_pvalue)] <- 0
write.table(wide_mean,'Th2_means.txt',quote=FALSE,sep='\t',row.names=FALSE)
write.table(wide_pvalue,'Th2_pvalue.txt',quote=FALSE,sep='\t',row.names=FALSE)
comm <- paste(mypython, ' plot dot_plot --means-path Th2_means.txt --pvalues-path Th2_pvalue.txt --output-name Th2_interaction_dotplot_by_cellphonedb.pdf --output-path ./', sep='')
system(comm)
###
metadata <- data.frame(CellType=data$CellType_b)
pdf('Th2_interaction_dotplot_by_ktplots_scale_mean.pdf',width=10,height=10)
p <- ktplots::plot_cpdb(
        scdata=metadata,
        cell_type1="Th2",
        cell_type2=".",  # this means all cell-types
        celltype_key="CellType",
        means=wide_mean,
        pvals=wide_pvalue,
        standard_scale=TRUE,
        max_size = 5
    )
print(p)
dev.off()


# Proliferating_T
data <- read.table('Proliferating_T.txt',header=TRUE,sep='\t',quote='')
data$CellType <- paste0(data$CellType_a,'|',data$CellType_b)
data$partner_a <- 'NULL'
data$partner_b <- 'NULL'
data$secreted <- 'NULL'
data$is_integrin <- 'NULL'
data$receptor_a <- 'NULL'
data$receptor_b <- 'NULL'
wide_mean <- dcast(data, id_cp_interaction+interacting_pair+partner_a+partner_b+gene_a+gene_b+secreted+receptor_a+receptor_b+annotation_strategy+is_integrin~ CellType, value.var = "mean")
wide_mean[is.na(wide_mean)] <- 0
wide_pvalue <- dcast(data, id_cp_interaction+interacting_pair+partner_a+partner_b+gene_a+gene_b+secreted+receptor_a+receptor_b+annotation_strategy+is_integrin~ CellType, value.var = "pValue")
wide_pvalue[is.na(wide_pvalue)] <- 0
write.table(wide_mean,'Proliferating_T_means.txt',quote=FALSE,sep='\t',row.names=FALSE)
write.table(wide_pvalue,'Proliferating_T_pvalue.txt',quote=FALSE,sep='\t',row.names=FALSE)
comm <- paste(mypython, ' plot dot_plot --means-path Proliferating_T_means.txt --pvalues-path Proliferating_T_means.txt --output-name Proliferating_T_interaction_dotplot_by_cellphonedb.pdf --output-path ./', sep='')
system(comm)
#
metadata <- data.frame(CellType=data$CellType_b)
pdf('Proliferating_T_interaction_dotplot_by_ktplots_scale_mean.pdf',width=10,height=10)
p <- ktplots::plot_cpdb(
        scdata=metadata,
        cell_type1="Proliferating T",
        cell_type2=".",  # this means all cell-types
        celltype_key="CellType",
        means=wide_mean,
        pvals=wide_pvalue,
        standard_scale=TRUE,
        max_size = 5
    )
print(p)
dev.off()


