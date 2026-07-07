## 输入一个sce对象和tissue名，得到top10 marker基因热图，保存top20基因csv

run_marker = function(sce,tissue){
  # sce是sce对象，可以是sce_list[[tissue]]
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
  DoHeatmap(sce.Scale,top_10)+ggtitle(tissue)
  
  ggsave(filename = paste0('./',tissue,'_top10_heatmap.png'),width = 12,height = 20)
  top20_df = apply(marker$names,2,head,20)
  write.csv(top20_df,file = paste0('./',tissue,'_top20_marker.csv'))
}