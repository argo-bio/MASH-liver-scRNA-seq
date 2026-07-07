## run cellchat
## input a sce obj and give out a cellchat obj

run_cellchat <- function(sce){
  ## creat cellchat object
  data.input <- sce@assays$RNA@layers$data
  colnames(data.input) = colnames(sce)
  rownames(data.input) = rownames(sce)
  meta <- sce@meta.data
  
  cellchat <- createCellChat(object = data.input,meta = meta, group.by = 'cell_type2')
  
  ## 载入数据库
  CellChatDB <- CellChatDB.mouse # 人类的就改成.human
  
  #CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中
  
  ## 分析计算
  #表达量预处理
  cellchat <- subsetData(cellchat,features = NULL)#取出表达数据
  cellchat <- identifyOverExpressedGenes(cellchat,do.fast = FALSE)#寻找高表达的基因
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj = PPI.mouse)#投影到PPI
  cellchat <- computeCommunProb(cellchat, raw.use = T)#默认计算方式为type = "truncatedMean",
  
  #默认cutoff的值为20%，即表达比例在25%以下的基因会被认为是0， trim = 0.1可以调整比例阈值
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  #去掉通讯数量很少的细胞
 
  cellchat <- computeCommunProbPathway(cellchat)
  #每对配受体的预测结果存在net中，每条通路的预测结果存在netp中
  cellchat <- aggregateNet(cellchat)
  #计算联路数与通讯概率，可用sources.use and targets.use指定来源与去向
}
