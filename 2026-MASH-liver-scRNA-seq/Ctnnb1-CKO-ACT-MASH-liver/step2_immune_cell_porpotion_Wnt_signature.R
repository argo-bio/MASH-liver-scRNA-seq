## June 11 2026

### step2 immune cell porpotion and wnt/beta-catenin signature score
calculate_cell_proportions <- function(sce) {
  
  meta_data <- as.data.frame(sce@meta.data)
  
  prop_data <- meta_data %>%
    group_by(group, cell_type_fine) %>%
    summarise(n_cells = n(), .groups = 'drop') %>%
    group_by(group) %>%
    mutate(total_cells = sum(n_cells),
           proportion = n_cells / total_cells)
  
  return(prop_data)
}

prop_data <- calculate_cell_proportions(sce)

ggplot(prop_data, aes(x = interaction(group), y = proportion, fill = cell_type_fine)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_brewer(palette = "Set3") +
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
ggsave('./output/part1/immune_cell_proportion.png',width = 7,height = 6,bg = 'white')
write.csv(prop_data,file = './output/part1/immune_cell_proportion.csv')


selCols = colorRampPalette(RColorBrewer::brewer.pal(n = 11, name = "Paired"))(11)
library(scRNAtoolVis)
p <- clusterCornerAxes(object = sce, reduction = 'umap',
                       clusterCol = "cell_type_fine", # 分组依据的列名
                       pSize =0.1,
                       arrowType ='open', # 坐标轴箭头类型
                       lineTextcol ='black',# 角线和标签颜色
                       cornerTextSize =5,
                       keySize =3, # 图注大小
                       show.legend =T, # 隐藏冗余图例
                       cellLabel =F,# 显示细胞标签
                       cellLabelSize =5,
                       noSplit =TRUE,
                       addCircle =F,# 开启亚群底纹云朵圈
                       cicAlpha =0.1,# 圆圈透明度
                       cicDelta =0.5,
                       nbin =200)
ggsave(paste0("./output/Fig4/immune_cell_umap.png"), plot = p, width = 8, height = 5)


wnt_gene = c("Apc", "Axin1", "Btrc", "Ccnd1", "Crebbp", "Csnk2a1", "Ctbp1", 
             "Ctnnb1", "Dvl1", "Frat1", "Fzd1", "Gsk3b", "Hdac1", "Hnf1a", 
             "Map3k7", "Myc", "Nlk", "Ppard", "Ppp2ca", "Smad4", "Tab1", 
             "Tle1")
wnt_gene %in% rownames(sce) %>% table()
gene=as.data.frame(rownames(sce))
wnt_gene = as.data.frame(wnt_gene)
wnt_gene = as.list(wnt_gene)

sce <- AddModuleScore(object = sce,features = wnt_gene,name="Wnt_Score")
VlnPlot(subset(sce, 
               subset = sce$cell_type_fine %in% macrophages), 
        features = "Wnt_Score1", group.by = "group",pt.size = 0.1)+geom_boxplot(width=.2,col="black",fill="white")
ggsave('./output/part1/all_mac_wnt_score.png',width = 6,height = 9)



sce$cell_group <- paste(sce$cell_type_fine, sce$group, sep = "_")
VlnPlot(subset(sce, 
               subset = sce$cell_type_fine %in% c( 'LAM','KC','mo-KC','KC-like LAM')), 
        features = "Wnt_Score1", group.by = "cell_group",pt.size = 0.1)+geom_boxplot(width=.2,col="black",fill="white")
ggsave('./output/part1/sub_mac_wnt_score.png',width = 18,height = 6)

VlnPlot(sce,features = "Wnt_Score1", group.by = "cell_type_fine",pt.size = 0.1,split.by = 'group')+
  geom_boxplot(width=.2,col="black",fill="white")