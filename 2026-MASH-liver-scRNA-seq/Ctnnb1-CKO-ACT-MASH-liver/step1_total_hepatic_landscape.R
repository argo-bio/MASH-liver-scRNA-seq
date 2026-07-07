# WT, CKO, ACT mice liver scRNA-seq

## Mar 18 2026

### step1 cell type annotation

library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(COSG)
library(ggplot2)

rm(list = ls())
sce <- readRDS('./rds/DZOE2025123793-b1_data_ob_v3.rds')

sce <- UpdateSeuratObject(sce)

DimPlot(sce,group.by = 'leiden')
DimPlot(sce,group.by = 'new_celltype')
DimPlot(sce,group.by = 'group')
Idents(sce) <- sce$leiden

marker = cosg(
  sce,
  groups='all',
  assay='RNA',
  slot='data',
  mu=1,
  n_genes_user=100,
  remove_lowly_expressed=TRUE,
  expressed_pct=0.1)#cosg快速计算单细胞marker

top_10 <- unique(as.character(apply(marker$names,2,head,10)))
sce.Scale <- ScaleData(sce ,features =  top_10)
DoHeatmap(sce.Scale,top_10)
ggsave('./output/leiden_top10_heatmap.png',width = 12,height = 30)

top20_df = apply(marker$names,2,head,20)
write.csv(top20_df,file = paste0('./output/leiden_top20_marker.csv'))

marker_list <- list(
  "Myeloid_General" = c("Lyz2", "Itgam", "Ccr2", "Cx3cr1","Ly6c2"),
  "Kupffer_Macrophage" = c("Adgre1", "Clec4f", "Timd4"),
  "LAM" = c("Trem2","Spp1","Fabp5","Marco","Folr2"),
  "T_Cells" = c("Cd3e", "Cd3g", "Cd4", "Cd8a", "Gzmb", "Foxp3"),
  "B_Cells" = c("Cd79a", "Cd19", "Ms4a1"),
  "cDC1" = c("Xcr1", "Itgax", "Clec9a"),
  "cDC2" = c("Cd209a","Mgl2"),
  "Neutrophils" = c("Ly6g", "S100a8", "S100a9", "Mmp9"),
  "NK_Cells" = c("Ncr1", "Klrb1c"),
  "pDC" = c("Siglech", "Ccr9", "Bst2"),
  "Endothelial" = c("Ptprc","Ptprb", "Kdr", "Flt4"),
  "Fibroblast" = c("Col1a1","Col3a1","Dcn"),
  "LESC" = c("Stab2", "Lyve1", "Fcgr2b"),
  "Hepatocytes" = c("Alb", "Hnf4a", "Cyp3a11")
)
dot_plot_simple <- DotPlot(sce, features = marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('./output/core_gene_expre.png',bg = 'white',width = 16,height = 6)

DimPlot(sce,group.by = 'leiden',label = T)

marker_list <- list(
  "Hep" = c("Alb", "Apoc3", "Pck1"),
  "Fibroblast" = c("Col1a1", "Dcn", "Pdgfra"),
  "Epi" = c("Epcam", "Cdh1", "Krt18"),
  "Basophils" = c("Cyp11a1","Ms4a2")
)
dot_plot_simple <- DotPlot(sce, features = marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('./output/non_immune_gene_expre.png',bg = 'white',width = 8,height = 6)

sce$cell_type <- case_when(
  sce$leiden %in% c(0,2,3,4,13,14,16,17) ~ "LESCs",
  sce$leiden %in% c(12) ~ "Hep",
  sce$leiden %in% c(22) ~ "Epi",
  sce$leiden %in% c(18) ~ "Fibroblast",
  sce$leiden %in% c(5,7,23,25) ~ "T_cells",
  sce$leiden %in% c(1,15,21,27) ~ "B_cells",
  sce$leiden %in% c(6,9,10,11,20,24) ~ "Mono-Macrophage",
  sce$leiden %in% c(8,26,19) ~ "Neutrophil",
  TRUE ~ "Unknown" # 处理未匹配的cluster
)

DimPlot(sce,group.by = 'cell_type',label = T)

# subset all immune cells
sce_immune <- subset(sce,cell_type %in% c("T_cells","B_cells","Mono-Macrophage","Neutrophil"))

sce_immune$group %>% table()
# saveRDS(sce_immune,file = './rds/all_immune_cells.rds')
