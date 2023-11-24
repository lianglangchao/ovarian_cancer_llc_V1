
library(Seurat)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggpubr)
library(pheatmap)
library(monocle)
library(RColorBrewer)

CD8T <- readRDS("CD8T_0209_anno_monocle3_sample30000.rds")

mtcell <- which(is.na(CD8T$percent.mt) == TRUE)

cell2 <- CD8T[,-mtcell]
data <- as(as.matrix(cell2@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = cell2@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
monocle_cds <- newCellDataSet(data,
                              phenoData = pd,
                              featureData = fd,
                              lowerDetectionLimit = 0.1,
                              expressionFamily = negbinomial.size())

HSMM <- monocle_cds
 <-M <- mateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
#Remi <- detectGenes(HSMM, min_expr = 3 )
print(head(fData(HSMM)))
print(head(pData(HSMM)))



##essed_genes <- cell2@assays$RNA@var.features

# Th_test_res <- differentialGeneTest(HSMM[expressed_genes,],
                                      fullModelFormulaStr = "~percent.mt")
ordeMordering_genes <- row.names(subset(diff_test_res, qval < 0.01))M
 <- setOrderingFilter(HSMM, ordering_genes)
plot_ordering_genes(HSMM)


#Tr <- reduceDimension(HSMM, max_components = 2,
   method = 'DDRTree')  <- orderCells(HSMM)

a <- pData(HSMM)
write.csv(a,"pdata.csv",quote=F)

pdf("CD8T_Pseudotime_0209.pdf",width = 6,height = 4)
plot_cell_trajectory(HSMM, color_by = "Pseudotime")
dev.off()


pdf("CD8T_celltype1_tissue_0209.pdf",width = 20,height = 4)
plot_cell_trajectory(HSMM, color_by = "celltype1")+ facet_wrap("~Tissue1", nrow = 1)
dev.off()

pdf("CD8T_celltype1_celltype1_0209.pdf",width = 20,height = 4)
plot_cell_trajectory(HSMM, color_by = "celltype1")+ facet_wrap("~celltype1", nrow = 1)
dev.off()

pdf("CD8T_state_0209.pdf",width = 6,height = 4)
plot_cell_trajectory(HSMM, color_by = "State")
dev.off()

saveRDS(HSMM,"CD8T_celltype1_monocle3_sample30000_monocle2_0209.rds")