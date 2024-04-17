setwd("d:/project/OV2/OV2_revision/")

library(Seurat)
endo2 <- readRDS("D:/project/OV2/09_endothelial_231127/OV-Endo-1127.rds")

library(RColorBrewer)
display.brewer.all()

cols1<-brewer.pal(12, "Set3")
mycolors1 <- colorRampPalette(cols1)(12)

Idents(endo2) <- "Tissue1"
endo3 <- subset(endo2,idents=c("Fallopian tube","Normal Ovarium","Ovarium"))

pdf("Fig3_1.pdf",width = 15,height = 6)
DimPlot(endo3,group.by = "celltype2",split.by = "Tissue1",cols = mycolors1,label = T)
dev.off()
