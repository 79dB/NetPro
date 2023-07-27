###5.disease 

library(data.table)
library(igraph)
library(openxlsx)
library(gplots)
library(pheatmap)
library(vegan)
library(ggplot2)
library(ggh4x)
library(tidyverse)

##泛癌部分
##药物靶点信息
setwd("github\\NetPro\\Data\\")
dat <- as.data.frame(fread("Drug target information.csv"))
##药物信息
drugs <- unique(dat$DrugBankID)
drugs <- sort(drugs)

##网络信息
PPI <- as.data.frame(fread("human_interactome_symbol.txt",header = F))
PPI <- PPI[,-3]
#画PPI网络图
g <- graph_from_data_frame(PPI,directed = F)
#判断图g是否联通
is_connected(g)
#rest为剩余节点
rest <- V(g)[components(g)$membership==1]
#将rest从igraph.vs数据类型转换为向量
rest <- as_ids(rest)
nodes <- unique(rest)

##筛选
PPI <- PPI[PPI$V1%in%nodes&PPI$V2%in%nodes,]
#画PPI网络图
g <- graph_from_data_frame(PPI,directed = F)
#判断图g是否联通
is_connected(g)

temp <- c()
for (i in 1:length(drugs)) {
  temp <- c(temp,nrow(dat[dat$DrugBankID%in%drugs[i],]))
}
##查看对应基因为1的药物,5个
min(temp)
length(which(temp==1))
which(temp==1)

##创建药物对应基因列表
drug_gene <- list()
for (i in 1:length(drugs)) {
  drug_gene[[i]] <- dat$Symbol[dat$DrugBankID%in%drugs[i]]
}
names(drug_gene) <- drugs

##读入疾病相关信息
disease <- as.data.frame(fread("pan-cancer.csv"))
##查看泛癌基因250个，交集242个
length(unique(disease$Gene))
length(intersect(unique(disease$Gene),nodes))
##保留最大子图中的基因
genes <- intersect(unique(disease$Gene),nodes)

setwd("github\\NetPro\\Data\\")
##读入Jaccard结果
result <- as.data.frame(fread("Side-reaction scores for strong synergy combinations.csv"))
result <- result[,-3]

##计算
dCC <- mean(shortest.paths(g,genes,genes))
for (i in 1:nrow(result)) {
  result$dAA[i] <- mean(shortest.paths(g,
                                       drug_gene[[result$drug_A[i]]],
                                       drug_gene[[result$drug_A[i]]]))
  result$dBB[i] <- mean(shortest.paths(g,
                                       drug_gene[[result$drug_B[i]]],
                                       drug_gene[[result$drug_B[i]]]))
  print(i)
}
result$dCC <- dCC

##计算dAC和dBC
for (i in 1:nrow(result)) {
  result$dAC[i] <- mean(shortest.paths(g,
                                       drug_gene[[result$drug_A[i]]],
                                       genes))
  result$dBC[i] <- mean(shortest.paths(g,
                                       drug_gene[[result$drug_B[i]]],
                                       genes))
  print(i)
}
##计算dAB
for (i in 1:nrow(result)) {
  result$dAB[i] <- mean(shortest.paths(g,
                                       drug_gene[[result$drug_A[i]]],
                                       drug_gene[[result$drug_B[i]]]))
  print(i)
}

##计算sAB,sAC和sBC以及score
for (i in 1:nrow(result)) {
  result$sAB[i] <- result$dAB[i] - (result$dAA[i] + result$dBB[i])/2
  result$sAC[i] <- result$dAC[i] - (result$dAA[i] + result$dCC[i])/2
  result$sBC[i] <- result$dBC[i] - (result$dBB[i] + result$dCC[i])/2
  print(i)
}
for (i in 1:nrow(result)) {
  result$score[i] <- (result$sAC[i]+result$sBC[i])/2
}

setwd("github\\NetPro\\Data\\")
##写出
#write.table(result,"Pan-cancer calculation results.csv",sep=",",col.names = T,row.names = F)

##保留三位小数进行处理，筛选前10%数据进行画图
##保留三位小数
for (i in 3:12) {
  result[,i] <- round(result[,i],3)
}
result10 <- result[,c(1,2,12)]
#result$score <- round(result$score,3)
result10 <- result[result$score<quantile(result$score,0.1),]
write.table(result10[,c(1,2,12)],"top 10% of pan-cancer result.csv",
            sep=",",col.names = T,row.names = F)

######
#2.单一药物部分
drug_info <- result10[,c(1,2,12)]
drug_info <- drug_info[,-3]
##查看药物数目,41
length(union(drug_info$drug_A,drug_info$drug_B))
##提取药物
drugs <- union(drug_info$drug_A,drug_info$drug_B)
drugs <- sort(drugs)

setwd("github\\NetPro\\Data\\")
##读入药物靶点信息
target_info <- as.data.frame(fread("Drug target information.csv"))
##产看药物数目,136
length(unique(target_info$DrugBankID))
##交集,41
length(intersect(unique(target_info$DrugBankID),union(drug_info$drug_A,drug_info$drug_B)))
##筛选对应药物的靶点信息
target_info <- target_info[target_info$DrugBankID%in%drugs,]
##剩余41个药物
length(unique(target_info$DrugBankID))

##构建药物靶点对应关系
drug_gene <- list()
for (i in 1:length(drugs)) {
  drug_gene[[i]] <- unique(target_info$Symbol[target_info$DrugBankID%in%drugs[i]])
}
names(drug_gene) <- drugs

##读入PPI网络
PPI <- as.data.frame(fread("human_interactome_symbol.txt",header = F))
PPI <- PPI[,-3]
#画PPI网络图
g <- graph_from_data_frame(PPI,directed = F)
#判断图g是否联通
is_connected(g)
#rest为剩余节点
rest <- V(g)[components(g)$membership==1]
#将rest从igraph.vs数据类型转换为向量
rest <- as_ids(rest)
nodes <- unique(rest)
##保留最大联通子图
PPI <- PPI[PPI$V1%in%nodes&PPI$V2%in%nodes,]
#画PPI网络图
g <- graph_from_data_frame(PPI,directed = F)
#判断图g是否联通
is_connected(g)

##筛选药物基因
for (i in 1:length(drug_gene)) {
  drug_gene[[i]] <- intersect(drug_gene[[i]],nodes)
}

##读入疾病相关信息
disease_info <- as.data.frame(fread("Disease location and corresponding gene.csv"))
##筛选10个以上基因的数据
disease_info <- disease_info[disease_info$num>=10,]
##去除神经和头颈信息
disease_info <- disease_info[disease_info$position!="nerve"&
                               disease_info$position!="head and neck",]
##疾病部位
diseases <- disease_info$position
disease_gene <- list()
for (i in 1:length(diseases)) {
  disease_gene[[i]] <- intersect(nodes,unique(unlist(strsplit(disease_info$gene[disease_info$position%in%diseases[i]],";"))))
}
names(disease_gene) <- diseases

##计算不同癌症的dAB
result <- list()
for (i in 1:length(diseases)) {
  result[[i]] <- drug_info
}
names(result) <- diseases

##计算dAA,dBB,dCC
dAA <- c()
for (i in 1:nrow(drug_info)) {
  dAA <- c(dAA,mean(shortest.paths(g,drug_gene[[drug_info$drug_A[i]]],
                                   drug_gene[[drug_info$drug_A[i]]])))
  print(i)
}
dBB <- c()
for (i in 1:nrow(drug_info)) {
  dBB <- c(dBB,mean(shortest.paths(g,drug_gene[[drug_info$drug_B[i]]],
                                   drug_gene[[drug_info$drug_B[i]]])))
  print(i)
}
dCC <- c()
for (i in 1:length(diseases)) {
  dCC <- c(dCC,mean(shortest.paths(g,disease_gene[[diseases[i]]],
                                   disease_gene[[diseases[i]]])))
  print(i)
}
##dAB
dAB <- c()
for (i in 1:nrow(drug_info)) {
  dAB <- c(dAB,mean(shortest.paths(g,drug_gene[[drug_info$drug_A[i]]],
                                   drug_gene[[drug_info$drug_B[i]]])))
  print(i)
}

##添加dAA,dBB,dCC
for (i in 1:length(result)) {
  result[[i]]$dAA <- dAA
  result[[i]]$dBB <- dBB
  result[[i]]$dCC <- dCC[i]
  result[[i]]$dAB <- dAB
}


##计算dAC,dBC
for (i in 1:length(result)) {
  for (j in 1:nrow(result[[i]])) {
    result[[i]]$dAC[j] <- mean(shortest.paths(g,drug_gene[[result[[i]]$drug_A[j]]],
                                              disease_gene[[diseases[i]]]))
    result[[i]]$dBC[j] <- mean(shortest.paths(g,drug_gene[[result[[i]]$drug_B[j]]],
                                              disease_gene[[diseases[i]]]))
    print(paste(c("i=",i),collapse = ""))
    print(paste(c("j=",j),collapse = ""))
  }
}
##计算sAB,sAC和sBC
for (i in 1:length(result)) {
  for (j in 1:nrow(result[[i]])) {
    result[[i]]$sAB[j] <- result[[i]]$dAB[j] - (result[[i]]$dAA[j]+
                                                  result[[i]]$dBB[j])/2
    result[[i]]$sAC[j] <- result[[i]]$dAC[j] - (result[[i]]$dAA[j]+
                                                  result[[i]]$dCC[j])/2
    result[[i]]$sBC[j] <- result[[i]]$dBC[j] - (result[[i]]$dBB[j]+
                                                  result[[i]]$dCC[j])/2
    print(paste(c("i=",i),collapse = ""))
    print(paste(c("j=",j),collapse = ""))
  }
}
##计算score
for (i in 1:length(result)) {
  for (j in 1:nrow(result[[i]])) {
    result[[i]]$score[j] <- (result[[i]]$sAC[j]+result[[i]]$sBC[j])/2
  }
}
bf <- result

##构建矩阵
output <- drug_info
##添加疾病得分的计算结果
for (i in 1:length(result)) {
  output <- cbind(output,result[[i]]$score)
}
colnames(output) <- c("drug1","drug2",diseases)

##写出
setwd("github\\NetPro\\Data\\")
write.table(output,"Matrix for single disease.csv",sep=",",col.names=T,row.names=F)

###3.热图部分
mat_info <- output
setwd("github\\NetPro\\Data\\")
##读入drugbank数据，获取药物名称
dbk <- as.data.frame(fread("drug links.csv"))
##转换名称
for (i in 1:nrow(mat_info)) {
  mat_info$drug1[i] <- dbk$Name[dbk$`DrugBank ID`%in%mat_info$drug1[i]]
  mat_info$drug2[i] <- dbk$Name[dbk$`DrugBank ID`%in%mat_info$drug2[i]]
}
##处理药物组合和行名
com <- c()
for (i in 1:nrow(mat_info)) {
  com <- c(com,paste(mat_info[i,1:2],collapse = "-"))
}
rownames(mat_info) <- com
mat_info <- mat_info[,-c(1,2)]
##列名大写
colnames(mat_info) <- paste(toupper(substr(colnames(mat_info),1,1)),
                            substr(colnames(mat_info),2,10),sep="")


##pheatmap热图
pheatmap(mat_info)
##进行z-score列转换
dat <- as.data.frame(scale(mat_info))
#出图
pheatmap(dat) 
##美化

# 使用clustering_method参数来设定不同聚类方法
#可以设定为'ward.D', 'ward.D2', 'single', 'complete', 'average', 
#'mcquitty', 'median' or 'centroid'
pheatmap(dat, clustering_method = "single")

# 行列聚类距离度量
#clustering_distance_rows="correlation"参数设定行聚类距离方法为Pearson corralation，
#默认为欧氏距离"euclidean"
pheatmap(dat,clustering_method = "median",
         clustering_distance_rows = "correlation",
         clustering_distance_cols= "correlation" )

# 去掉边框线
# 去掉边框线
pheatmap(dat,clustering_method = "centroid",
         clustering_distance_rows = "correlation",
         clustering_distance_cols= "correlation",
         border=FALSE)

# 数值小于-1的用*表示
##保存热图 10*13
pheatmap(dat, clustering_method = "single",
         border=FALSE,angle_col = 45,
         display_numbers = matrix(ifelse(dat < -1, "*"," "), nrow(dat)),
         fontsize_number = 15,
         fontsize_col = 13,fontsize_row = 13)

