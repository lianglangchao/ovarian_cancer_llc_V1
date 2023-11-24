
library(Seurat)
library(ggplot2)
library(RColorBrewer)

TC <- readRDS("T_0205_recluster.rds")

TC1 <- TC[,which(TC@assays$RNA@counts["CD3D",] >0)]
TC1 <- NormalizeData(TC1)
TC1 <- FindVariableFeatures(TC1, selection.method = "vst", nfeatures = 2000)
TC1 <- ScaleData(TC1, verbose = T)
TC1 <- RunPCA(TC1,features = VariableFeatures(object = TC1))


library(harmony)
TC1 <- RunHarmony(TC1, group.by.vars = "Patient_source",max.iter.harmony=15)
TC1 <- RunUMAP(TC1, reduction = "harmony", dims = 1:30)
TC1 <- FindNeighbors(TC1, reduction = "harmony", dims = 1:30) %>% FindClusters()

saveRDS(TC1,"T_0208_CD3D_recluster.rds")
