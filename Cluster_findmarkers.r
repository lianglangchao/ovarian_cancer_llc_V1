
library(Seurat)
library(ggplot2)


myOC <- readRDS("OvC-harmony-0201_filtered.rds")

table(myOC$seurat_clusters)

Idents(myOC) <- myOC$major
major.markers <- FindAllMarkers(myOC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(major.markers,"major_ct_allmarkers.csv",quote = F)

Idents(myOC) <- myOC$RNA_snn_res.0.8
oc.markers <- FindAllMarkers(myOC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(oc.markers,"cluster_allmarkers.csv",quote = F)
