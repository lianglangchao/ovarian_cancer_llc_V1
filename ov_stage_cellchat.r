
library(CellChat)
library(patchwork)
library(Seurat)
library(ggpubr)
library(cowplot)
options(stringsAsFactors = FALSE)


####创建CellChat 对象

ov1 <- readRDS("OvC-harmony-0201_filtered.rds")
ov1

Idents(ov1) <- "Tissue1"
ov2 <- subset(ov1,idents="Ovarium")
ov2
remove(ov1)
gc()

ov2$Stage2 <- ov2$Stage1
ov2$Stage2[which(ov2$Stage1 %in% c("IA","IC1","IC2"))] = "Stage-I"
ov2$Stage2[which(ov2$Stage1 %in% c("IIB"))] = "Stage-II"
ov2$Stage2[which(ov2$Stage1 %in% c("IIIB","IIIC"))] = "Stage-III"
ov2$Stage2[which(ov2$Stage1 %in% c("IV","IVA","IVB"))] = "Stage-IV"
ov2$Stage2[which(is.na(ov2$Stage1))] = "No information"

Idents(ov2) <- "Stage2"
ov3 <- subset(ov2,idents=c("Cancer-free","Stage-I","Stage-II","Stage-III","Stage-IV"))

remove(ov2)
gc()

saveRDS(ov3,"ov_stage2.rds")


for (i in unique(ov3$Stage2)) {
  
  cell <- subset(ov3,idents=i)
  cell$major <- as.character(cell$major)
  data.input  <- cell
  meta <- cell@meta.data
  unique(meta$major)
  
  cellchat <- createCellChat(object = data.input, meta = meta,group.by = "major")
  cellchat
  CellChatDB <- CellChatDB.human 
  
  
  CellChatDB.use <- CellChatDB 
  cellchat@DB <- CellChatDB.use
  
  
  cellchat <- subsetData(cellchat) 
  future::plan("multiprocess", workers = 5) # do parallel
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.human)
  
  
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  table(cellchat@idents)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  saveRDS(cellchat, paste("ov",i,"cellchat.rds",sep = "_"))  
}


library(CellChat)
library(patchwork)
library(Seurat)
library(ggpubr)
library(cowplot)
options(stringsAsFactors = FALSE)


data.dir <- './comparison'
dir.create(data.dir)
setwd(data.dir)

cellchatC <- readRDS("../ov_Cancer-free_cellchat.rds")
cellchatI <- readRDS("../ov_Stage-I_cellchat.rds")
cellchatII <- readRDS("../ov_Stage-II_cellchat.rds")
cellchatIII <- readRDS("../ov_Stage-III_cellchat.rds")
cellchatIV <- readRDS("../ov_Stage-IV_cellchat.rds")

object.list <- list(Cancerfree = cellchatC, StageI = cellchatI,
                    StageII = cellchatII,StageIII = cellchatIII,
                    StageIV = cellchatIV)
cellchat <- mergeCellChat(object.list, add.names = names(object.list),cell.prefix = T)
cellchat

pdf("compareInteractions.pdf",width = 10,height = 5)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1:5))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1:5), measure = "weight")
gg1 + gg2
dev.off()


net1 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(3,8),  comparison = c(1,2,3,4,5), angle.x = 45)
write.csv(net1$data,"comparison_source_bubble_data.csv",quote = F)

pdf("comparison_source_bubble.pdf",width = 6,height = 35)
net1
dev.off()

net2 <- netVisual_bubble(cellchat, sources.use = 4, targets.use = c(3,8),  comparison = c(1,2,3,4,5), angle.x = 45,signaling = "MHC-I")

pdf("comparison_source_bubble_MHC-I.pdf",width = 8,height = 5)
net2
dev.off()


interaction_name2 <- net2$data$interaction_name[-which(net2$data$interaction_name_2 %in% c("HLA-F - CD8B","HLA-F - CD8A"))]
interaction_name2 <- as.data.frame(interaction_name2)

pdf("comparison_source_bubble_MHC-I.pdf",width = 8,height = 4)
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(3,8),  comparison = c(1,2,3,4,5), angle.x = 45,pairLR.use = interaction_name2)
dev.off()


net1 <- netVisual_bubble(cellchat, sources.use = 8, targets.use = c(3,4),  comparison = c(1,2,3,4,5), angle.x = 45)
write.csv(net1$data,"T_comparison_source_bubble_data.csv",quote = F)

pdf("T_comparison_source_bubble.pdf",width = 6,height = 7)
net1
dev.off()


net1 <- netVisual_bubble(cellchat, sources.use = 3, targets.use = c(4,8),  comparison = c(1,2,3,4,5), angle.x = 45)
write.csv(net1$data,"Epi_comparison_source_bubble_data.csv",quote = F)

pdf("Epi_comparison_source_bubble.pdf",width = 6,height = 30)
net1
dev.off()

net1 <- netVisual_bubble(cellchat, sources.use = c(1:8), targets.use = 8,  comparison = c(1,2,3,4,5), angle.x = 45)
write.csv(net1$data,"All_T_comparison_source_bubble_data.csv",quote = F)

pdf("All_T_comparison_source_bubble.pdf",width = 12,height = 30)
net1
dev.off()

pdf("All_T_comparison_source_bubble_TIGIT.pdf",width = 12,height = 3)
netVisual_bubble(cellchat, sources.use = c(2,3,4,6,7), targets.use = 8,  comparison = c(1,2,3,4,5), angle.x = 45,signaling = c("NECTIN","PVR"))
dev.off()

