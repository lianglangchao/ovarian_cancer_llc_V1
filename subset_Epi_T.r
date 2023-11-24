
library(Seurat)
library(ggplot2)
library(RColorBrewer)

myOC <- readRDS("OvC-harmony-0201_filtered.rds")
table(myOC$seurat_clusters)

Idents(myOC) <- "major"
unique(myOC$major)
cell <- subset(myOC,idents="Epithelial cell")
saveRDS(cell,"Epithelial_0205.rds")

cell <- subset(myOC,idents="T cell")
saveRDS(cell,"T_0205.rds")

cell <- subset(myOC,idents="Fibroblast")
saveRDS(cell,"Fibroblast_0205.rds")