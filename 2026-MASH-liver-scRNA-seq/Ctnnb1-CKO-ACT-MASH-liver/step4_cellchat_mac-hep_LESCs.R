## Jun 25 2026

### step4 CellChat between macrophages and hepatocytes and LESCs
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(COSG)
library(ggplot2)
library(stringr)
library(ggpubr)
library(ggimage)
library(SingleCellExperiment)

rm(list = ls())

### merge marcophage and LESC/Hepatocytes
sce <- readRDS('./rds/immune_cell_fine_annotation.rds')
DimPlot(sce,group.by = 'cell_type_fine',label = T)

imm <- sce
rm(sce)

hep <- subset(sce,sce$leiden %in% c(0,12))
DimPlot(hep,label = T)

imm$cell_type_fine %>% unique()
hep$cell_type_fine <- ifelse(hep$leiden==0,'LESC','Hepatocytes')

hep@meta.data=hep@meta.data[,c(1:3,8,28)]
imm@meta.data=imm@meta.data[,c(1:3,5,8)]

sce <- merge(hep,imm)

ct <- CreateSeuratObject(counts = cbind(sce@assays$RNA@layers$counts.1,sce@assays$RNA@layers$counts.2),
                         meta.data = sce@meta.data)
colnames(ct) <- colnames(sce)
rownames(ct) <- rownames(sce)
rm(hep,imm,sce)

colname <- colnames(sce)
rowname <- rownames(sce)

sce <- ct
rm(ct)

source('./script/run_harmony.R')
sce$orig.ident <- sce$group
sce <- run_harmony(sce,0.8)

DimPlot(sce,group.by = 'seurat_clusters',label = T)
DimPlot(sce,group.by = 'group',label = T)
DimPlot(sce,group.by = 'cell_type_fine',label = T)
saveRDS(sce,file = './rds/cellchat.rds')

### merged cell clusters UMAP
sce <- readRDS('./rds/cellchat.rds')
df <- sce@meta.data[,c(5,7)]
df$cell <- paste0(df$cell_type_fine,'-',df$seurat_clusters)

sce$cell_seu <- paste0(sce$cell_type_fine,'-',sce$seurat_clusters)
ne_df <- table(sce$cell_seu) %>% as.data.frame()
ne_df <- ne_df[ne_df$Freq>100,]

sce$cell_type <- case_when(
  sce$seurat_clusters %in% c(0,25,7,13) ~ "B_cell",
  sce$seurat_clusters %in% c(16) ~ "cDC1",
  sce$seurat_clusters %in% c(10) ~ "Hepatocytes",
  sce$seurat_clusters %in% c(3,9,20,14) ~ "KC",
  sce$seurat_clusters %in% c(18,21) ~ "LAM",
  sce$seurat_clusters %in% c(1,2) ~ "LESC",
  sce$seurat_clusters %in% c(11,12,15) ~ "MDM",
  sce$seurat_clusters %in% c(5) ~ "Monocyte",
  sce$seurat_clusters %in% c(6) ~ "NK",
  sce$seurat_clusters %in% c(4,8) ~ "T_cell",
  sce$seurat_clusters %in% c(17) ~ "pDC",
  sce$seurat_clusters %in% c(19) ~ "Mast",
  sce$seurat_clusters %in% c(22,24) ~ "Neutrophil",
  TRUE ~ "Unknown" # 处理未匹配的cluster
)

DimPlot(sce,group.by = 'cell_type',label = T)

rownames(sce) <- rowname
colnames(sce) <- colname
sce <- subset(sce,sce$cell_type!='Unknown')

sce$cell_type2 <- ifelse(sce$cell_type %in% c('KC','LAM','MDM'),'macrophage',sce$cell_type)
DimPlot(sce,group.by = 'cell_type2',label = T)

library(RColorBrewer)
my_cols <- c(brewer.pal(11, "Set3"), c("#526E2D99","#FD744699"))

sce$cell_type <- factor(sce$cell_type,levels = c('B_cell',"cDC1","KC","LAM","Mast","MDM","Monocyte","Neutrophil","NK","pDC","T_cell","Hepatocytes","LESC"))
DimPlot(sce,group.by = 'cell_type',label = F) + 
  scale_color_manual(values = my_cols)
ggsave('./output/cellchat/UAMP.pdf',width = 8,height = 6)
ggsave('./output/cellchat/UAMP.png',width = 8,height = 6)

### cell chat
library(CellChat)
library(cowplot)
library(patchwork)

sce.act <- subset(sce,sce$group=='ACT')
sce.wt <- subset(sce,sce$group=='WT')
sce.cko <- subset(sce,sce$group=='CKO')

source('./script/run_cellchat.R')

cellchat.wt <- run_cellchat(sce.wt)

netVisual_bubble(cellchat.wt,sources.use = 5,targets.use = 3:4)

cellchat.act <- run_cellchat(sce.act)
cellchat.cko <- run_cellchat(sce.cko)

save(cellchat.wt, cellchat.act, cellchat.cko,
     file = "./rds/cell_chat_output.RData")
load("./rds/cell_chat_output.RData")
object.list <- list(WT=cellchat.wt,ACT=cellchat.act,CKO=cellchat.cko)
cellchat <- mergeCellChat(object.list,add.names = names(object.list))

netVisual_bubble(cellchat,sources.use = 2:3,targets.use = 5,comparison = c(1,2,3))

df <- CellChatDB[["interaction"]][["interaction_name"]] %>% as.data.frame()
df2 <- CellChatDB[["interaction"]][["interaction_name_2"]] %>% as.data.frame()

ggsave('./output/cellchat/cellchat.pdf',width = 8,height = 10)

netVisual_bubble(cellchat,sources.use = 5,targets.use = 1:4,comparison = c(1,2,3))
ggsave('./output/cellchat/cellchat_mac_to_hep.pdf',width = 8,height = 10)
