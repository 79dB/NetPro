#######github项目
##4.Side reaction

library(data.table)
library(randomForest)
library(varSelRF)
library(pROC)
library(ROCR)
library(caret)
library(modEvA)
library(e1071)
library(nnet)
library(openxlsx)

##1.预测协同和拮抗组合
setwd("github\\NetPro\\Data\\")
##读入金标准
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
seed <- 431
set.seed(seed)
neg <- neg[sample(nrow(neg),nrow(pos)),]
set <- rbind(pos,neg)

##读入所有候选组合
etwd("github\\NetPro\\Data\\")
cand <- as.data.frame(fread("The filtered network propagation score.csv"))
##寻找非金标准部分
com_std <- c()
for (i in 1:nrow(std)) {
  com_std <- c(com_std,paste(std[i,1:2],collapse = ";"))
}
com_cand <- c()
for (i in 1:nrow(cand)) {
  com_cand <- c(com_cand,paste(cand[i,1:2],collapse = ";"))
}
length(intersect(com_cand,com_std))
#位置
cand <- cand[-which(com_cand%in%com_std),]
##设置种子
seed <- 431
set.seed(seed)
##训练集和测试集
train <- set
test <- cand
##寻找最优的mtry参数,最优为2
rate=1     #设置模型误判率向量初始值

for(i in 1:8){
  set.seed(seed)
  rf_train<-randomForest(class~.,data=train[,-c(1,2)],mtry=i,ntree=1000)
  rate[i] <- mean(rf_train$err.rate)   #计算基于OOB数据的模型误判率均值
  #print(rf_train)    
}

rate 
plot(rate)

##模型训练
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,-c(1,2)],mtry=2,
                         ntree=1000,
                         important=TRUE,proximity=TRUE)

##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,-c(1,2)])#"prob"显示概率
rf_prob <- predict(rf_model,newdata=test[,-c(1,2)],"prob")
table(rf_pred)

##添加类别
cand$class <- rf_pred
##添加概率
cand <- cbind(cand,rf_prob)
##区分协同和拮抗组合
pos <- cand[cand$class=="synergistic",]
pos <- pos[,-c(11,12)]
neg <- cand[cand$class=="antagonistic",]
neg <- neg[,-c(11,13)]


##添加协同类别
for (i in 1:nrow(pos)) {
  if(pos$synergistic[i]>=0.5&pos$synergistic[i]<0.6){
    pos$grade[i] <- "Slight"
  }else if(pos$synergistic[i]>=0.6&pos$synergistic[i]<0.8){
    pos$grade[i] <- "Moderate"
  }else
    pos$grade[i] <- "Strong"
}

##添加拮抗类别
for (i in 1:nrow(neg)) {
  if(neg$antagonistic[i]>=0.5&neg$antagonistic[i]<0.6){
    neg$grade[i] <- "Slight"
  }else if(neg$antagonistic[i]>=0.6&neg$antagonistic[i]<0.8){
    neg$grade[i] <- "Moderate"
  }else
    neg$grade[i] <- "Strong"
}
##更改列名
colnames(pos) <- c(colnames(pos)[1:10],"Synergistic probability","Synergistic strength")
colnames(neg) <- c(colnames(pos)[1:10],"Antagonistic probability","Antagonistic strength")
colnames(pos) <- c("drug_A","drug_B",colnames(pos)[3:12])
colnames(neg) <- c("drug_A","drug_B",colnames(neg)[3:12])

##2.副反应得分计算
setwd("github\\NetPro\\Data\\")
dat <- as.data.frame(fread("ADReCS.csv"))
##筛选强协同
xt <- pos[pos$`Synergistic probability`>=0.8,]

length(unique(dat$drugbank))    ##1362
length(union(xt$drug_A,xt$drug_B))    ####强协同剩余607  

##筛选严重类别为轻微的数据,剩余17277
dat <- dat[dat$class%in%"Mild",]
##剩余药物,579
length(unique(dat$drugbank))
##药物
d1 <- unique(dat$drugbank)
d2 <- union(xt$drug_A,xt$drug_B)
##强协同剩余251   
drug_ins <- intersect(d1,d2)
drug_ins <- sort(drug_ins)

##筛选副反应信息中包含对应药物的数据,强协同剩余7487  
dat <- dat[dat$drugbank%in%drug_ins,]
##筛选协同信息中包含对应药物的数据,强协同剩余735 
xt <- xt[xt$drug_A%in%drug_ins&xt$drug_B%in%drug_ins,]
##4强协同剩余151   
length(union(xt$drug_A,xt$drug_B))

##构建协同药物对应副反应列表
adr <- list()
drug <- union(xt$drug_A,xt$drug_B)
drug <- sort(drug)
for (i in 1:length(drug)) {
  adr[[i]] <- unique(dat$ADR_TERM[dat$drugbank%in%drug[i]])
}
names(adr) <- drug
##统计副反应
n <- c()
for (i in 1:length(adr)) {
  n <- c(n,length(adr[[i]]))
}
##强协同副反应最多191
max(n)
##所有不良反应
adr_all <- c()
for (i in 1:length(adr)) {
  adr_all <- c(adr_all,adr[[i]])
}
adr_all <- unique(adr_all)    ##强协同剩余2172ADR 
##统计药物对应副反应
a <- data.frame(drug=drug)
a$freq <- n

##计算Jaccard系数
result <- xt[,c(1,2)]
for (i in 1:nrow(result)) {
  result$Jaccard[i] <- length(intersect(adr[[result$drug_A[i]]],adr[[result$drug_B[i]]]))/length(union(adr[[result$drug_A[i]]],adr[[result$drug_B[i]]]))
}
##保留两位小数
result$Jaccard <- round(result$Jaccard,2)

##强协同平均值0.01
##强协同剩余575个组合,136个药物 
result <- result[result$Jaccard<=mean(result$Jaccard),]
length(union(result$drug_A,result$drug_B))

##写出
setwd("github\\NetPro\\Data\\")
write.table(result,"Side-reaction scores for strong synergy combinations.csv",
            sep=",",col.names = T,row.names = F)
