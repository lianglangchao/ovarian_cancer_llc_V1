
library(Seurat)
library(ggplot2)
library(harmony)
library(RColorBrewer)

EPI <- readRDS("EPI_0206_recluster.rds")

EPI <- FindClusters(EPI,resolution = 0.3)
EPI$'RNA_snn_res.0.3' <- factor(EPI$'RNA_snn_res.0.3',levels = 0:19)

all_marker <- FindAllMarkers(EPI,only.pos = T,logfc.threshold = 0.25,min.diff.pct = 0.25)
write.csv(all_marker,"All_EPI_clu_DEGs.csv")