
threshold <- 0


my_vln <- function(sce,gene,sub_celltype){
  p <- VlnPlot(subset(sce, subset = GetAssayData(sce, assay = "RNA")[gene, ] > 0 & 
                   sce$cell_type %in% sub_celltype), 
          features = gene, group.by = "new_group")+geom_boxplot(width=.2,col="black",fill="white")&
    stat_compare_means(method="t.test",hide.ns = F,
                       comparisons = c(list(c("sd", "wd-12w")),
                                       list(c("sd", "wd-24w")),
                                       list(c("sd", "wd-36w"))),
                       label="p.signif", 
                       bracket.size=0.8,      
                       tip.length=0,        
                       size=6)&  
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
  
  expr=GetAssayData(sce, assay="RNA", layer='data')
  expr_df <- expr[gene,]
  plot_data <- data.frame(
    cell_id = colnames(sce),
    expression = expr_df,
    group = sce$new_group,
    celltype = sce$cell_type
  ) %>%
    mutate(
      # 判断是否为阳性细胞
      positive = expression > threshold,
      # 创建分组+细胞类型的组合变量
      group_celltype = paste(group, celltype, sep = "_")
    )
  plot_data = plot_data[plot_data$celltype %in% sub_celltype,]
  positive_stats <- plot_data %>%
    group_by(group) %>%
    summarise(
      total_cells = n(),
      positive_cells = sum(positive),
      positive_ratio = positive_cells / total_cells
    )
  return(list(plot = p, stats = positive_stats))
}



