

library(Seurat)
library(ggplot2)
library(RColorBrewer)

TC1 <- readRDS("Tcell_0209_anno.rds")

color1 <- colorRampPalette(c(brewer.pal(8, "Accent"),brewer.pal(8, "Set2")))(20)

da1 <- TC1@meta.data
da1 <- da1[da1$Tissue1 == "Ovarium",]
da2 <- as.data.frame(table(da1$celltype1,da1$Stage1))
colnames(da2)[1:2] <- c("Cluster","Source")
p <- ggplot(da2,mapping = aes(Source,Freq,fill=Cluster))+
  geom_bar(stat="identity",width = 0.6,position="fill")+
  scale_fill_manual(values =  color1[1:7])+
  theme_bw()+
  #coord_flip()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),      #去除边界线
        axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle=0,vjust = 0,hjust = 0.5))#加上坐标线       

pdf("Ovarium_stage1_celltype1_barplot.pdf",width = 7,height = 3)
print(p)
dev.off()

