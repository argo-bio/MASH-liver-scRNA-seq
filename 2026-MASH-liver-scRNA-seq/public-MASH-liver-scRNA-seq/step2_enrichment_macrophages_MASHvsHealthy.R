## Aug 25 2025

### step2 porpotion and enrichment

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
library(tidyr)
library(scales)
rm(list = ls())

sce <- readRDS('./rds/all_immune_cells.rds')

sce$fine_group <- ifelse(sce$group == "sd", "sd",
                         paste0("wd","-", sce$stage))
sce$fine_group <- ifelse(sce$group == "wd", sce$fine_group,
                         paste0("sd","-", sce$stage))
sce$fine_group <- factor(sce$fine_group, 
                         levels = c("sd-12w", "sd-24w","sd-36w","wd-12w", "wd-24w", "wd-36w"))
sce$sd_group <- ifelse(sce$group == "sd", "sd",
                       paste0("wd","-", sce$stage))
sce$sd_group <- factor(sce$sd_group, 
                       levels = c("sd", "wd-12w", "wd-24w", "wd-36w"))

## porpotion
calculate_cell_proportions <- function(sce) {
  meta_data <- as.data.frame(sce@meta.data)
  prop_data <- meta_data %>%
    group_by(fine_group, cell_type) %>%
    summarise(n_cells = n(), .groups = 'drop') %>%
    group_by(fine_group) %>%
    mutate(total_cells = sum(n_cells),
           proportion = n_cells / total_cells)
  return(prop_data)
}

prop_data <- calculate_cell_proportions(sce)
df <- read.csv('./output/recluster_myeloid/immune_cell_proportion.csv',row.names = 'X')
prop_data$proportion <- df$proportion
ggplot(prop_data, aes(x = interaction(fine_group), y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(  values = c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", 
                                 "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", 
                                 "#CCEBC5", "#FFED6F")) +
  labs(title = "",
       x = "",
       y = "Proportion",
       fill = "Cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(size = 14))
ggsave('./output/recluster_myeloid/immune_cell_proportion.png',width = 7,height = 6,bg = 'white')
#write.csv(prop_data,file = './output/recluster_myeloid/immune_cell_proportion.csv')

# activate Wnt Signature

wnt_activation_genes <- c("Ctnnb1", "Ccnd1", "Myc", "Axin2","Lef1","Tcf4","Fzd1","Fzd2","Dvl2")
gene=as.data.frame(rownames(sce))
wnt_gene = as.data.frame(wnt_activation_genes)
wnt_gene = as.list(wnt_gene)
sce <- AddModuleScore(object = sce,features = wnt_gene,name="Wnt_Score_act")

wnt_repress_genes <- c("Apc","Axin1", "Gsk3b","Btrc","Ckb1")
wnt_gene = as.data.frame(wnt_repress_genes)
wnt_gene = as.list(wnt_gene)
sce <- AddModuleScore(object = sce,features = wnt_gene,name="Wnt_Score_rep")

sce$Wnt_Score1 <- sce$Wnt_Score_act1 - sce$Wnt_Score_rep1

sce$cell_type_group <- paste0(sce$cell_type,'-',sce$group)
VlnPlot(subset(sce, subset = sce$group %in% c('wd') & sce$cell_type %in% c('MDM','KC','LAM')), 
        features = "Wnt_Score1", group.by = "cell_type",pt.size = 0.1)+geom_boxplot(width=.2,col="black",fill="white")&
  stat_compare_means(method="t.test",hide.ns = F,
                     comparisons = c(list(c("MDM", "LAM"),c("MDM","KC"),c("KC",'LAM'))),
                     label="p.signif", 
                     bracket.size=0.8,      
                     tip.length=0,        
                     size=6)&  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave('./output/recluster_myeloid/all_mac_celltype_wnt_score.png',width = 6,height = 9)

# Enrichment analyze
library(Seurat)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(msigdbr)
library(clusterProfiler)



# Visulization
selCols = colorRampPalette(RColorBrewer::brewer.pal(n = 12, name = "Paired"))(12)

# devtools::install_github("junjunlab/scRNAtoolVis")
library(scRNAtoolVis)
p <- clusterCornerAxes(object = sce, reduction = 'umap',
                       clusterCol = "cell_type",
                       pSize =0.1,
                       arrowType ='open', 
                       lineTextcol ='black',
                       cornerTextSize =5,
                       keySize =1.5, 
                       show.legend =T, 
                       cellLabel =F,
                       cellLabelSize =5,
                       noSplit =TRUE,
                       addCircle =F,
                       cicAlpha =0.1,
                       cicDelta =0.5,
                       nbin =200) +
  scale_color_manual(values = selCols)
ggsave(paste0("./output/Fig1/immune_cell_umap.png"), plot = p, width = 8, height = 5)


ggplot(prop_data, aes(x = interaction(fine_group), y = proportion, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(  values = selCols) +
  labs(title = "",
       x = "",
       y = "Proportion",
       fill = "Cell type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        legend.position = "right",
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        axis.text.y = element_text(size = 14))
ggsave(paste0("./output/Fig1/immune_cell_proportion.pdf"), width = 6, height = 5)
ggsave(paste0("./output/Fig1/immune_cell_proportion.png"), width = 6, height = 5,bg = 'white')


VlnPlot(subset(sce, subset = sce$cell_type %in% macrophage), 
        features = "Wnt_Score_act1", group.by = "group",pt.size = 0)+
  geom_boxplot(width=.2,col="black",fill="white")&
  stat_compare_means(method="t.test",hide.ns = F,
                     comparisons = c(list(c("sd", "wd"))),
                     label="p.signif", 
                     bracket.size=0.8,      
                     tip.length=0,        
                     size=10)&  
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
ggsave(paste0("./output/Fig1/macropahge_wnt_score.png"), width = 6, height = 5,bg = 'white')


marker_list <- list(
  "T_Cells" = c("Cd3e", "Cd3g", "Cd8a"),
  "NK_Cells" = c("Ncr1", "Klrb1c"),
  "B_Cells" = c("Cd79a", "Cd19", "Ms4a1"),
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
ggsave('./output/Fig1/all_immune_marker_dotplot.pdf',bg = 'white',width = 14,height = 6)

