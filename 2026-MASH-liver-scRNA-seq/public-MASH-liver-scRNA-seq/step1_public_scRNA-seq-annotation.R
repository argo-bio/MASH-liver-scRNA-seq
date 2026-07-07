# public avaliable MASH mice liver scRNA-seq (GSE156059)

## Aug 25 2025

### step1 cell type annotation

library(dplyr)
library(Seurat)
library(patchwork)

# arrange files into folder
dir <- './data/GSE208750/'
fs <- list.files(dir,'genes.tsv.gz')
samples1 <- gsub('_genes.tsv.gz','',fs)
library(stringr)
samples2=str_split(samples1,'_',simplify = T)[,1]
setwd(dir)
lapply(1:length(samples2), function(i){
  x=samples2[i]
  y=samples1[i]
  dir.create(x,recursive = T) #sample ID as fold name
  file.copy(from=paste0(y,'_genes.tsv.gz'), #file name
            to=file.path(x,  'genes.tsv.gz' )) 
  file.copy(from=paste0(y,'_matrix.mtx.gz'),
            to= file.path(x, 'matrix.mtx.gz' ) ) 
  file.copy(from=paste0(y,'_barcodes.tsv.gz'),
            to= file.path(x, 'barcodes.tsv.gz' )) 
})


setwd('../GSE208750-Liver/')
files <- list.files('./',full.names = T)
sceList = lapply(files, function(patient){
  print(patient)
  ct <- Read10X(patient)
  sce <- CreateSeuratObject(counts = ct,
                            project = str_split(patient,pattern = '/')[[1]][2],
                            min.cells=3,
                            min.features = 200)
  return(sce)
})
sce <- merge(x=sceList[[1]],y=sceList[-1])
sce <- JoinLayers(sce)
table(sce$orig.ident)

sce <- JoinLayers(sce)
saveRDS(sce,file = './rds/sce_all_merge_readin.rds')

### annotation
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(COSG)
library(ggplot2)

rm(list = ls())
sce <- readRDS('./rds/sce_all_merge_readin.rds')

source('./script/run_harmony.R')
sce <- run_harmony(sce)

sce$seurat_clusters %>% table()
sce@meta.data %>% colnames()

DimPlot(sce,group.by = 'seurat_clusters',reduction = 'umap',label = T)
ggsave('./output/step2cluster_marker/dim_plot.png',width = 8,height = 6)

source('./script/run_marker.R')
run_marker(sce,'MAFLD_liver')

if (T) {
  Endo = c('Pecam1','Cdh5','Cldn5')
  Fibroblasts = c('Col1a1', 'Col3a1', 'Col5a1', 'Dcn', 'Twist1')
  Epi = c('Epcam','Krt8','Vil1')
  ProliferatingCells = c('Mki67','Birc5','Top2a')
  
  Immune = c('PTPRC')
  
  BCells = c('Cd79a','Cd19','Ms4a1')
  Plasma = c('Jchain','Mzb1','Sdc1')
  
  Mac = c('Lyz2','Cd68','Cx3cr1')
  Monocyte = c('Ly6i', 'Ly6c2', 'Ccr2', 'Ms4a4c', 'Ass1')
  Mast = c('Cma1','Hdc','Tpsb2')
  Neutrophils = c('S100a8','S100a9','Csf3r')
  pDCs = c('Siglech','Irf7','Ccr9')
  cDCs=c('Xcr1', 'Flt3',  'Ccr7')
  
  TandNK = c('Cd3e','Cd3d','Cd3g')
  
  CD8T = c('Cd8b1', 'Cd8a', 'Gzmk', 'Themis', 'Ifng', 'Cxcr6')
  CD4Treg = c('Foxp3', 'Ctla4', 'Il2ra', 'Tnfrsf18')
  CD4Th = c('Cd3g', 'Cd3d', 'Trac', 'Cd3e', 'Cd4', 'Icos', 'Pik3cd', 'Bcl2', 'Bcl11b', 'Itk')
  
  cDC1 = c("CLEC9A","XCR1","CLNK","CADM1","ENPP1","SNX22","NCALD","DBN1","HLA-DOB","PPY")
  cDC2= c( "CD1C","FCER1A","CD1E","AL138899.1","CD2","GPAT3","CCND2","ENHO","PKIB","CD1B")
  cDC3 = c("HMSD","ANKRD33B","LAD1","CCR7","LAMP3","CCL19","CCL22","INSM1","TNNT2","TUBB2B")
  
  M2_Mac = c('Mrc1', 'Csf1r', 'Trem2', 'Maf', 'Ccl12', 'Rnase4', 'Gatm', 'Ckb', 'Cd72')
  M1_Mac = c('Nos2', 'Ptges', 'Il1a', 'F10', 'Ptgs2', 'SIc7a2', 'Il7r', 'Il1b')
  
  
  marker_list1 = list(
    Endo = Endo,
    Epi = Epi,
    Fibroblasts = Fibroblasts,
    Prolife = ProliferatingCells,
    Immu = Immune,
    BCell = BCells,
    Plasma = Plasma,
    Mac = Mac,
    Mono = Monocyte,
    Mast = Mast,
    Neutro = Neutrophils,
    pDCs = pDCs,
    cDCs = cDCs,
    CD8T = CD8T,
    CD4Treg = CD4Treg,
    CD4Th = CD4Th
  ) #classical markers 
}

marker_list1 <- lapply(marker_list1, function(x){
  marker_list1 = str_to_title(x)
  return(marker_list1)
})

DotPlot(sce, features = marker_list1 )  + 
  theme(axis.text.x=element_text(angle=45,hjust = 1))+ggtitle('MAFLD_liver')

ggsave('./marker_dot_plot.png',bg = 'white',width = 20,height = 8)

sce$cell_type <- case_when(
  sce$seurat_clusters %in% c(8,14) ~ "pDC",
  sce$seurat_clusters %in% c(6,11) ~ "cDC2",
  sce$seurat_clusters %in% c(10,19) ~ "cDC1",
  sce$seurat_clusters %in% c(4,18) ~ "Hep",
  sce$seurat_clusters %in% c(15) ~ "Epi",
  sce$seurat_clusters %in% c(2) ~ "Monocyte",
  sce$seurat_clusters %in% c(3,13) ~ "Myelocyte",
  sce$seurat_clusters %in% c(7) ~ "MDM",
  sce$seurat_clusters %in% c(20) ~ "Neutrophil",
  sce$seurat_clusters %in% c(1,12,17) ~ "KC",
  sce$seurat_clusters %in% c(0) ~ "LAM",
  sce$seurat_clusters %in% c(5,16) ~ "T_cell",
  sce$seurat_clusters %in% c(9) ~ "B_cell",
  TRUE ~ "Unknown" 
)

DimPlot(sce,group.by = 'cell_type',label = T)
ggsave('./output/recluster_myeloid/myeloid_umap.png',width = 8,height = 6)

marker_list <- list(
  "Monocyte" = c("Ly6c2","Ccr2","Cx3cr1",
                 "Itgam"),
  "KC" = c("Adgre1", "Clec4f", "Timd4"),
  "LAM" = c("Trem2","Spp1","Fabp5"),
  "pDC" = c("Siglech", "Ccr9", "Bst2"),
  "cDC1" = c("Xcr1","Clec9a","Itgax"),
  "cDC2" = c("Cd209a","Mgl2"),
  "Neutrophils" = c("Ly6g","S100a8", "S100a9", "Mmp9")
)
DotPlot(sce, features = marker_list,group.by = 'cell_type') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('./output/recluster_myeloid/celltype_marker.png',bg = 'white',width = 9,height = 6)



