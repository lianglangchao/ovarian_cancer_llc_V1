setwd("/dellfsqd2/ST_LBI/USER/chaichaochao/Graduation/02_data_plot/03_harmony_V2/02_annotation/03_Fibro/")

library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
FB <- readRDS("Fibro_0316_recluster.rds")

color <- colorRampPalette(c(brewer.pal(8, "Accent"),brewer.pal(8, "Set3"),brewer.pal(8, "Dark2")))(27)

da1 <- FB@meta.data
da3 <- da1[-which(da1$Source1 %in% c("OvC16","OvC20","OvC22")),]
da4 <- da3[which(da3$Tissue1 == "Ovarium"),]
FBOV <- FB[,rownames(da4)]

Idents(FBOV) <- "RNA_snn_res.0.3"
myoFB <- subset(FBOV,idents=0)

# ECMgene <- c("ACTA2","MYL9","POSTN","COL10A1","COL11A1","MMP11","COMP","TAGLN","FN1")


myoFB$Stage1[is.na(myoFB$Stage1)] <- "No information"

# pdf("myoFB_Ovarium_ECM_dotplot.pdf",height=5,width = 8)
# DotPlot(myoFB,features = ECMgene,group.by = "Stage1")+coord_flip()
# dev.off()

ECMgene1 <- read.table("GO0031012_ECM.txt",header = T,sep ="\t")
myoFB <- AddModuleScore(myoFB,features = list(ECMgene1$gene),name = "ECM")

pdf("myoFB_Ovarium_ECM_vlnplot_0918.pdf",width = 8,height = 4)
VlnPlot(myoFB,features = "ECM1",pt.size = 0,group.by = "Stage1")
dev.off()

data <- data.frame(Stage=myoFB$Stage1, score=myoFB$ECM1)
data <- data[-which(data$Stage == "No information"),]
library(RColorBrewer)
Stage_col <- c(rev(brewer.pal(6, "BuGn")[3:6]),brewer.pal(6, "YlOrRd"))
names(Stage_col) <- unique(data$Stage)[order(unique(data$Stage))]

compaired <- list(c("IC2", "IC1"),c("IIB", "IC1"),c("IIIB", "IC1"),c("IIIC", "IC1"),
                  c("IV", "IC1"),c("IVA", "IC1"),c("IVB", "IC1"))
p <- ggplot(data,aes(Stage,score,fill=Stage))+
  geom_violin(scale = "width")+
  scale_fill_manual(values = Stage_col)+
  theme(plot.title=element_text(size = 15),axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),panel.grid=element_blank(),legend.key = element_blank())+labs(x='Stage', y= 'score')+
  ggsignif::geom_signif(comparisons = compaired,step_increase = 0.1,
                        map_signif_level = T,test = t.test)+theme_bw()
pdf("myoFB_Ovarium_ECM_vlnplot_0918_1.pdf",width = 9,height = 5)
print(p)
dev.off()


p <- ggplot(data,aes(Stage,score,fill=Stage))+
  geom_violin(scale = "width")+
  scale_fill_manual(values = Stage_col)+
  theme(plot.title=element_text(size = 15),axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),axis.title.x=element_text(size = 15),
        axis.title.y=element_text(size = 15),panel.grid=element_blank(),legend.key = element_blank())+labs(x='Stage', y= 'score')+
  ggsignif::geom_signif(comparisons = compaired,step_increase = 0.1,
                        map_signif_level = F,test = t.test)+theme_bw()
pdf("myoFB_Ovarium_ECM_vlnplot_0918_2.pdf",width = 9,height = 5)
print(p)
dev.off()