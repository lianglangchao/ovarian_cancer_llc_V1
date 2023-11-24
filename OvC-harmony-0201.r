
library(Seurat)
library(ggplot2)
library(harmony)

ov1 <- readRDS("OvC-merge-0107.rds")
ov1$Patient_source[which(ov1$'Patient-source' %in% unique(ov1$'Patient-source')[-58])] <- ov1$'Patient-source'
unique(ov1$Patient_source)

ov1 <- ov1[,which(ov1$nFeature_RNA > 500)]

ov1 <- RunHarmony(ov1, group.by.vars = "Patient_source")
ov1 <- RunUMAP(ov1, reduction = "harmony", dims = 1:30)
ov1 <- FindNeighbors(ov1, reduction = "harmony", dims = 1:30) %>% FindClusters()

saveRDS(ov1,"OvC-harmony-0201.rds")

