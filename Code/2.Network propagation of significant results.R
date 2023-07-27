#######github项目
##2.Network propagation

library(data.table)
library(igraph)
library(RandomWalkRestartMH)
library(ggplot2)

##读入带有PCA的测度得分
setwd("github\\NetPro\\Data\\")
PPI <- as.data.frame(fread("Drug measure with PCA score.csv"))
##保留PCA得分
PPI <- PPI[,-c(9,10)]

##计算z-score和显著性p值
PPI$Z_score <- (PPI$PCA - mean(PPI$PCA))/sd(PPI$PCA)
PPI$p_value <- 2*pnorm(PPI$Z_score,lower.tail = F)
##筛选显著的结果
PPI <- PPI[PPI$p_value<0.05,]
##校正p值
PPI$p_adj <- p.adjust(PPI$p_value,method = "fdr")
##不叫正剩余45298，矫正后剩余17376

##将所有药物作为种子
seeds <- union(unique(PPI$Drug1),unique(PPI$Drug2))
##将种子进行排序
seeds <- sort(seeds)

##结果网络图和性质分析
specific <- PPI
##生成网络
PPI_Network <- graph.data.frame(specific,directed=FALSE)
##边权重
E(PPI_Network)$weight <- specific$PCA

##网络简化，删除多边和环
PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)

#RWR
PPI_MultiplexObject <- create.multiplex(list(PPI_Network),Layers_Name=c("PPI"))
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)
##模块挖掘
module <- list()
for (i in 1:length(seeds)) {
  module[[i]] <- Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,PPI_MultiplexObject,seeds[i])[[1]]
  print(i)
}
names(module) <- seeds

##添加种子
for (i in 1:length(module)) {
  module[[i]] <- cbind(seed=rep(seeds[i],nrow(module[[i]])),module[[i]])
}
##备份
bf <- module
##筛选网络中可以与种子对应上的节点
for (i in 1:length(module)) {
  module[[i]] <- module[[i]][module[[i]]$NodeNames%in%PPI$Drug2[PPI$Drug1%in%names(module)[i]],]
}

##去除空列表
n <- c()
for (i in 1:length(module)) {
  n <- c(n,nrow(module[[i]]))
}
sum(n)
which(n==0)
module[which(n==0)] <- NULL
##剩余1325个模块

##拼接结果
result <- data.frame()
for (i in 1:length(module)) {
  result <- rbind(result,module[[i]])
  print(i)
}
##查看药物数目,1440
length(union(result$seed,result$NodeNames))

rwr <- result
colnames(rwr) <- c(colnames(result)[-3],"RWR")
##筛选大于中位数,0.002
rwr <- rwr[rwr$RWR>median(rwr$RWR),]

##添加其他信息
num <- c()
for (i in 1:nrow(rwr)) {
  num <- c(num,which(PPI$Drug1%in%rwr$seed[i]&PPI$Drug2%in%rwr$NodeNames[i]))
  print(i)
}
rwr$ATC <- PPI$ATC[num]
rwr$Structure <- PPI$Structure[num]
rwr$GO <- PPI$GO[num]
rwr$KEGG <- PPI$KEGG[num]
rwr$Jaccard <- PPI$Jaccard[num]
rwr$Distance <- PPI$Distance[num]
rwr$PCA <- PPI$PCA[num]
rwr$z_score <- PPI$Z_score[num]
rwr$p_value <- PPI$p_value[num]
rwr$p_adj <- PPI$p_adj[num]

##调整顺序
rwr <- rwr[,c(1,2,5,4,9,8,6,7,10,3)]
write.table(rwr,"The filtered network propagation score.csv",
            sep=",",col.names = T,row.names = F)

