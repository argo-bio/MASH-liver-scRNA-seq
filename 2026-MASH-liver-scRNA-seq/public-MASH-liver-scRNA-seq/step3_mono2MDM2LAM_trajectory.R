## Apr 19 2026

### step3 poseto-trajectory of mono-macrophage

rm(list = ls())
sce <- readRDS('./rds/all_immune_cells.rds')
macro_all <- subset(sce, 
                    subset = cell_type %in% c("MDM", "KC", "LAM","Monocyte"))
DimPlot(macro_all,group.by = 'cell_type',label = T)
sce <- macro_all
sce <- CreateSeuratObject(counts = macro_all@assays$RNA@layers$counts,
                          project = 'MASH_mac',
                          min.cells=0,
                          min.features = 0)
colnames(sce)  <- colnames(macro_all)
rownames(sce) <- rownames(macro_all)
sce$group <- macro_all$group
sce$stage <- macro_all$stage
sce$cell_type_original <- macro_all$cell_type
sce$orig.ident <- sce$group

source('./script/run_harmony.R')
sce <- run_harmony(sce,0.8)

DimPlot(sce,label = T)
ggsave('./output/mac_pesodo/mac_umap.png',width = 8,height = 6)

marker_list <- list(
  "Mye" = c("Ptprc","Lyz2","Adgre1"),
  "Monocyte" = c("Ly6c2","Ccr2","Cx3cr1","Itgam"),
  "kc" = c("Cd5l", "Vsig4", "Slc1a3", "Cd163", "Gfra2", 
           "Adrb1", "Vcam1", "Hmox1","Timd4"),
  "lam" = c("Gpnmb", "Trem2", "Gdf15", "Ccl9", "Cd63", 
            "Fam20c", "Pla2g7", "Cd9", "Lgals3", "Cpp1", 
            "Fam20p5", "Lgals1", "Mmp12","Acp5", "Fbp1", "Spp1")
)

DotPlot(sce, features = marker_list,group.by = 'seurat_clusters') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('./output/mac_pesodo/mac_celltype_marker.png',bg = 'white',width = 14,height = 6)

DimPlot(sce,group.by = 'cell_type_original')

source('./script/run_marker.R')
run_marker(sce,'output/mac_pesodo/mac')
sce$cell_type <- case_when(
  sce$seurat_clusters %in% c(4,10,14,18) ~ "Monocyte",
  sce$seurat_clusters %in% c(1,2,11,5,7) ~ "MDM",
  sce$seurat_clusters %in% c(0,6,8,9,12,15,16) ~ "KC",
  sce$seurat_clusters %in% c(17) ~ "KC-like LAM",
  sce$seurat_clusters %in% c(3,13) ~ "Trem2+LAM",
  sce$seurat_clusters == 19 ~ "other",
  TRUE ~ "Unknown" # 处理未匹配的cluster
)
DimPlot(sce,group.by = 'cell_type')
saveRDS(sce,'./rds/mac_sub_celltype.rds')
sce <- readRDS('./rds/mac_sub_celltype.rds')
sce <- subset(sce,cell_type != 'other')

DimPlot(sce,group.by = 'cell_type')
ggsave('./output/mac_pesodo/mac_celltype_umap.png',width = 8,height = 6)
ggsave('./output/mac_pesodo/mac_celltype_uma.pdf',width = 8,height = 6)

marker_list <- list(
  "KC" = c("Lyz2","Adgre1","Vsig4","Timd4"),
  "Monocyte" = c("Ly6c2","Ccr2","Cx3cr1","Itgam"),
  "LAM" = c("Trem2","Gpnmb","Cd9","Spp1","Fabp5")
)
DotPlot(sce, features = marker_list,group.by = 'cell_type') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave('./output/mac_pesodo/mac_celltype_marker_validate.png',bg = 'white',width = 9,height = 5)
ggsave('./output/mac_pesodo/mac_celltype_marker_validate.pdf',bg = 'white',width = 9,height = 5)


# Monocle2
#devtools::install_local('../../monocle.tar.gz')
# remote install of monocle

library(monocle)
library(Seurat)
library(dplyr)
library(igraph)
packageVersion('monocle')
packageVersion('Seurat')
packageVersion('dplyr')

sce <- readRDS('../data/mac_sub_celltype.rds')

data <- as(as.matrix(sce@assays$RNA@layers$counts), 'sparseMatrix')
data[1:4,1:4]
rownames(data) <- rownames(sce)
colnames(data) <- colnames(sce)
pd <- new('AnnotatedDataFrame', data = sce@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=10, relative_expr = TRUE)

disp_table <- dispersionTable(mycds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
mycds <- setOrderingFilter(mycds, disp.genes)
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
suppressWarnings(mycds <- orderCells(mycds))

saveRDS(mycds,'../output/mono_mac_monocle.rds')

mycds <- readRDS('../output/mono_mac_monocle.rds')
mycds$cell_type = ifelse(mycds$cell_type=='other','Monocyte',mycds$cell_type)
mycds$group_cell_type <- paste0(mycds$group,'-',mycds$stage,'-',mycds$cell_type)
mycds$cell_type <- factor(mycds$cell_type,levels = c('KC','KC-like LAM','MDM','Monocyte','Trem2+LAM'))
plot_cell_trajectory(mycds, color_by = "cell_type") +
  facet_wrap(~group_cell_type, nrow = 6)
ggsave('../output/momo-mac-trajectory.png',width = 15,height = 25)
ggsave('../output/momo-mac-trajectory.pdf',width = 15,height = 25)
plot_cell_trajectory(mycds, color_by = "cell_type")|
  plot_cell_trajectory(mycds, color_by = "State")|
  plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave('../output/momo-mac-total-trajectory.png',width = 15,height = 6)
ggsave('../output/momo-mac-total-trajectory.pdf',width = 15,height = 6)

mycds$fine_group = paste0(mycds$group,mycds$stage)
plot_cell_trajectory(mycds, color_by = "cell_type")+
  facet_wrap(~fine_group, nrow = 2)

pData(mycds)$Trem2 = log2(exprs(mycds)["Trem2",] +1)
plot_cell_trajectory(mycds,color_by = "Trem2") +
  scale_color_continuous(type = "viridis")
ggsave('../output/momo-mac-Trem2-trajectory.png',width = 8,height = 6)
ggsave('../output/momo-mac-Trem2-trajectory.pdf',width = 8,height = 6)

# BEAM heatmap
BEAM_res <- BEAM(mycds, branch_point = 2, cores = 10, progenitor_method = 'duplicate')#查看分支点2两侧的基因表达变化
saveRDS(BEAM_res,'../output/BEAM_res.rds')
BEAM_res <- readRDS('../output/BEAM_res.rds')
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_geme <- row.names(subset(BEAM_res,qval < 0.01))
length(BEAM_geme)

library(ClusterGVis)

df <- plot_genes_branched_heatmap2(mycds[BEAM_geme2,],
                                   branch_point = 2,
                                   num_clusters = 4,
                                   cores = 1,
                                   use_gene_short_name = T,
                                   show_rownames = T)
visCluster(obj = df,plotType  = "heatmap",markGenes = marker_genes[marker_genes$cluster=='Trem2+LAM',]$gene)


df <- plot_genes_branched_heatmap2(mycds[BEAM_geme2,mycds$cell_type%in%c('Monocyte','Trem2+LAM','MDM')],
                                   branch_point = 2,
                                   num_clusters = 4,
                                   cores = 1,
                                   use_gene_short_name = T,
                                   show_rownames = T)
visCluster(obj = df,plotType  = "heatmap",markGenes = marker_genes[marker_genes$cluster=='MDM',]$gene)

BEAM_res <- BEAM_res[order(BEAM_res$qval),]

df <- plot_genes_branched_heatmap2(mycds[sig_genes,mycds$cell_type%in%c('Monocyte','Trem2+LAM','MDM')&mycds$Pseudotime<15.5],
                                   branch_point = 2,
                                   num_clusters = 4,
                                   cores = 1,
                                   use_gene_short_name = T,
                                   show_rownames = T)
gene <- c('Marco','Vsig4','Mrc1',
          'Gpnmb','Spp1','Ms4a7','Fabp5','Trem2','Cd9','Lamp2','Tgfbr1','Sirpa','Slpi','Slc11a1',
          'Rpl10','H3f3a','Rpl39','Rpl7a',
          'Ccr2','Chil3','Ly6c2','Clec4e','Fn1','S100a8','S100a9','Mki67')

pdf('./output/mono-mac-trajectory/branch-gene-heatmap2.pdf',height = 12,width = 9)
visCluster(obj = df,plotType  = "heatmap",markGenes = gene)
dev.off()