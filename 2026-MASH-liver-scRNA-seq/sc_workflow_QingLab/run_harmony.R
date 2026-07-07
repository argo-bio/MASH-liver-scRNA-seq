# 自定义函数：harmony整合、降维和聚类
## 输入一个sce对象，得到聚类后的结果

run_harmony <- function(input_sce,res){
  print(dim(input_sce))
  input_sce <- NormalizeData(input_sce, 
                             normalization.method = "LogNormalize",
                             scale.factor = 1e4) 
  input_sce <- FindVariableFeatures(input_sce)
  input_sce <- ScaleData(input_sce)
  input_sce <- RunPCA(input_sce, features = VariableFeatures(object = input_sce))
  seuratObj <- RunHarmony(input_sce, "orig.ident")
  names(seuratObj@reductions)
  seuratObj <- RunUMAP(seuratObj,  dims = 1:15, 
                       reduction = "harmony")
  input_sce=seuratObj
  input_sce <- FindNeighbors(input_sce, reduction = "harmony",
                             dims = 1:15) 
  input_sce.all=input_sce
  input_sce.all=FindClusters(input_sce.all, #graph.name = "CCA_snn", 
                             resolution = res, algorithm = 1)
  table(input_sce.all@active.ident)
  return(input_sce.all)
}