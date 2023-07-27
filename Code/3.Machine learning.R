#######github项目
##3.Machine learning

library(data.table)
library(randomForest)
library(varSelRF)
library(pROC)
library(ROCR)
library(caret)
library(modEvA)
library(e1071)
library(nnet)

########################################
#1.Calculate the efficiency scores of different classifiers
setwd("github\\NetPro\\Data\\")

##函数
performance<-function(table,n=2){
  if(!all(dim(table)==c(2,2)))
    stop("Must be a 2 X 2 table")
  tn=table[1,1]
  fp=table[1,2]
  fn=table[2,1]
  tp=table[2,2]
  sensitivity=tp/(tp+fn)
  ppv=tp/(tp+fp)
  accuracy=(tp+tn)/(tp+tn+fp+fn)
  result<-paste("Sensitivity = ",round(sensitivity,n),
                "\nPositive Predictive Value = ",round(ppv,n),
                "\nAccuracy = ",round(accuracy,n),"\n",seq=" ")
  cat(result)
}

##读入金标准
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)

##设置结果,顺序为：RF、SVM、Bys、Log、ANN、整合投票
result <- data.frame(classification=c("RF","SVM","Bys","Log","ANN","Integrated"))

##RF分类器相关信息
##设置种子
seed <- 431
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,-c(1,2)],mtry=2,
                         ntree=1000,
                         important=TRUE,proximity=TRUE)

##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,-c(1,2)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,-c(1,2)],"prob")
##生成四格表
rf_perf<-table(test$class,rf_pred,dnn=c("Actual","Predicted"))
rf_perf
test_perf <- rf_perf
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[1] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[1] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
##伪阳性率
result$FPR[1] <- 1-result$Specificity[1]
##假发现率
result$FDR[1] <- 1-result$Sensitivity[1]
##准确度
result$Accuracy[1] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
##阳性预测率
result$PPV[1] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
##阴性预测率
result$NPV[1] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
##F1评分
result$F1[1] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[1] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
##基尼指数
result$Gini[1] <- 2*result$Accuracy[1] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[1] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))


##SVM分类器相关信息
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
##设置种子
seed <- 322
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
svm_model <- svm(class~.,data = train[,-c(1,2)],kernel="radial",
                 cost=1,gamma=1/ncol(train[,-c(1,2)]),probability = TRUE)
##SVM预测
svm_pred <- predict(svm_model,newdata = test[,-c(1,2)])
svm_prob <- attr(predict(svm_model,newdata = test[,-c(1,2)],probability = T), "probabilities")
##生成四格表
svm_perf<-table(test$class,svm_pred,dnn=c("Actual","Predicted"))
svm_perf
test_perf <- svm_perf
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[2] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[2] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
##伪阳性率
result$FPR[2] <- 1-result$Specificity[2]
##假发现率
result$FDR[2] <- 1-result$Sensitivity[2]
##准确度
result$Accuracy[2] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
##阳性预测率
result$PPV[2] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
##阴性预测率
result$NPV[2] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
##F1评分
result$F1[2] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[2] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
##基尼指数
result$Gini[2] <- 2*result$Accuracy[2] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[2] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))


##Bys分类器相关信息
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)

seed <- 623
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
bys_model <- naiveBayes(class~. ,data = train[,-c(1,2)])
##贝叶斯预测
bys_pred <- predict(bys_model,newdata = test[,-c(1,2)])
##提取分类概率
bys_prob <- predict(bys_model,newdata = test[,-c(1,2)],type = "raw")
##生成四格表
bys_perf<-table(test$class,bys_pred,dnn=c("Actual","Predicted"))
bys_perf
test_perf <- bys_perf
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[3] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[3] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
##伪阳性率
result$FPR[3] <- 1-result$Specificity[3]
##假发现率
result$FDR[3] <- 1-result$Sensitivity[3]
##准确度
result$Accuracy[3] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
##阳性预测率
result$PPV[3] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
##阴性预测率
result$NPV[3] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
##F1评分
result$F1[3] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[3] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
##基尼指数
result$Gini[3] <- 2*result$Accuracy[3] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[3] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))


##Log分类器相关信息
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
seed <- 37
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
log_model <- glm(class~.,data = train[,-c(1,2)],family = "binomial")
log_prob <- predict(object =log_model,newdata=test[,-c(1,2)],type = "response")
log_pred<-ifelse(log_prob>=0.5,"synergistic","antagonistic")
log_pred <- factor(log_pred)
##四格表
log_perf <- table(test$class,log_pred,dnn=c("Actual","Predicted"))
log_perf
test_perf <- log_perf
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[4] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[4] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
##伪阳性率
result$FPR[4] <- 1-result$Specificity[4]
##假发现率
result$FDR[4] <- 1-result$Sensitivity[4]
##准确度
result$Accuracy[4] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
##阳性预测率
result$PPV[4] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
##阴性预测率
result$NPV[4] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
##F1评分
result$F1[4] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[4] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
##基尼指数
result$Gini[4] <- 2*result$Accuracy[4] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[4] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))


##ANN分类器相关信息
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
seed <- 150
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
##模型
ann_model <- nnet(class~.,data = train[,-c(1,2)],size=10)
#测试集
ann_pred <- predict(ann_model,newdata = test[,-c(1,2)],type="class")
ann_prob <- predict(ann_model,newdata = test[,-c(1,2)])
ann_perf <- table(test$class,ann_pred,dnn=c("Actual","Predicted"))
ann_perf
test_perf <- ann_perf
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[5] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[5] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
##伪阳性率
result$FPR[5] <- 1-result$Specificity[5]
##假发现率
result$FDR[5] <- 1-result$Sensitivity[5]
##准确度
result$Accuracy[5] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
##阳性预测率
result$PPV[5] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
##阴性预测率
result$NPV[5] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
##F1评分
result$F1[5] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[5] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
##基尼指数
result$Gini[5] <- 2*result$Accuracy[5] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[5] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))


##整合投票相关信息
setwd("github\\NetPro\\Data\\")
##设置种子
seed <- 431
set.seed(seed)
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
seed <- 431
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##数据
set <- rbind(pos,neg)

##SVM
seed <- 3
set.seed(seed)
svm_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
svm_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(svm_train)))),]
set.seed(seed)
svm_model <- svm(class~.,data = svm_train[,-c(1,2)],kernel="radial",
                 cost=1,gamma=1/ncol(svm_train[,-c(1,2)]),probability = TRUE)
svm_pred <- predict(svm_model,newdata = svm_test[,-c(1,2)])

##Logistic分类器
seed <- 3
set.seed(seed)
log_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
log_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(log_train)))),]
set.seed(seed)
log_model <- glm(class~.,data = log_train[,-c(1,2)],family = "binomial")
log_prob <- predict(object =log_model,newdata=log_test[,-c(1,2)],type = "response")
log_pred<-ifelse(log_prob>=0.5,"synergistic","antagonistic")
log_pred <- factor(log_pred)

##ANN分类器
seed <- 6
set.seed(seed)
ann_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
ann_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(ann_train)))),]
set.seed(seed)
ann_model <- nnet(class~.,data = ann_train[,-c(1,2)],size=10)
ann_pred <- predict(ann_model,newdata = ann_test[,-c(1,2)],type="class")
ann_pred <- factor(ann_pred)
names(ann_pred) <- rownames(ann_test)

##NB分类器
seed <- 6
set.seed(seed)
bys_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
bys_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(bys_train)))),]
set.seed(seed)
bys_model <- naiveBayes(class~. ,data = bys_train[,-c(1,2)])
bys_pred <- predict(bys_model,newdata = bys_test[,-c(1,2)])
names(bys_pred) <- rownames(bys_test)

##交集信息
test_ins <- c(rownames(svm_test),rownames(bys_test),
              rownames(log_test),rownames(ann_test))
test_ins <- as.data.frame(table(test_ins))
test_ins$test_ins <- as.character(test_ins$test_ins)
test_ins <- test_ins[test_ins$Freq>=3,]

##添加类别
for (i in 1:nrow(test_ins)) {
  if(length(svm_pred[names(svm_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$SVM[i] <- svm_pred[names(svm_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$SVM[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(bys_pred[names(bys_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$Bys[i] <- bys_pred[names(bys_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$Bys[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(log_pred[names(log_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$Log[i] <- log_pred[names(log_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$Log[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(ann_pred[names(ann_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$ANN[i] <- ann_pred[names(ann_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$ANN[i] <- NA
  
}
##去除频率信息
test_ins <- test_ins[,-2]
##投票判断
for (i in 1:nrow(test_ins)) {
  if(length(which(test_ins[i,-1]%in%2))>=3){
    test_ins$pred_class[i] <- "synergistic"
  }else
    test_ins$pred_class[i] <- "antagonistic"
}
##添加正式类别
for (i in 1:nrow(test_ins)) {
  test_ins$actual_class[i] <- as.character(set$class[which(rownames(set)%in%test_ins$test_ins[i])])
}
test_ins$pred_class <- factor(test_ins$pred_class)
test_ins$actual_class <- factor(test_ins$actual_class)
##四格表
int_perf <- table(test_ins$actual_class,test_ins$pred_class,
                  dnn=c("Actual","Predicted"))
int_perf
test_perf <- int_perf
##计算效能得分
##灵敏度,又叫真阳性率
result$Sensitivity[6] <- test_perf[2,2]/(test_perf[2,2]+test_perf[2,1])
##特异度,又叫真阴性率
result$Specificity[6] <- test_perf[1,1]/(test_perf[1,1]+test_perf[1,2])
##伪阳性率
result$FPR[6] <- 1-result$Specificity[6]
##假发现率
result$FDR[6] <- 1-result$Sensitivity[6]
##准确度
result$Accuracy[6] <- (test_perf[1,1]+test_perf[2,2])/sum(test_perf)
##阳性预测率
result$PPV[6] <- test_perf[2,2]/(test_perf[2,2]+test_perf[1,2])
##阴性预测率
result$NPV[6] <- test_perf[1,1]/(test_perf[1,1]+test_perf[2,1])
##F1评分
result$F1[6] <- (2*test_perf[2,2])/(2*test_perf[2,2]+test_perf[1,2]+test_perf[2,1])
##Mathews系数(MCC)，即Phi相关系数
result$MCC[6] <- (test_perf[2,2]*test_perf[1,1]-test_perf[1,2]*test_perf[2,1])/sqrt((test_perf[2,2]+test_perf[1,2])*(test_perf[2,2]+test_perf[2,1])*(test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1]))
##基尼指数
result$Gini[6] <- 2*result$Accuracy[6] - 1
##Kappa,通常0.75以上表示一致性较满意，0.4以下一致性不好
result$Kappa[6] <- ((test_perf[1,1]+test_perf[1,2])*(test_perf[1,1]+test_perf[2,1])+(test_perf[2,1]+test_perf[2,2])*(test_perf[1,2]+test_perf[2,2]))/(sum(test_perf)*sum(test_perf))

#########################################
#2.Draw AUC curves for different classifiers
##读入金标准
setwd("github\\NetPro\\Data\\")
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)

##设置种子
seed <- 431
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,-c(1,2)],mtry=2,
                         ntree=1000,
                         important=TRUE,proximity=TRUE)

##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,-c(1,2)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,-c(1,2)],"prob")

rf_roc <- roc(test$class,rf_prob[,2])
plot(rf_roc,print.auc=T)
plot(smooth(rf_roc),print.auc=T)

##SVM
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
##设置种子
seed <- 322
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]


set.seed(seed)
svm_model <- svm(class~.,data = train[,-c(1,2)],kernel="radial",
                 cost=1,gamma=1/ncol(train[,-c(1,2)]),probability = TRUE)
##SVM预测
svm_pred <- predict(svm_model,newdata = test[,-c(1,2)])
svm_prob <- attr(predict(svm_model,newdata = test[,-c(1,2)],probability = T), "probabilities")
##生成四格表
svm_perf<-table(test$class,svm_pred,dnn=c("Actual","Predicted"))
svm_perf
##查看
performance(svm_perf)
##ROC
svm_roc <- roc(test$class,svm_prob[,1])
plot(svm_roc,print.auc=T)
plot(smooth(svm_roc),print.auc=T)

##朴素贝叶斯
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)

seed <- 623
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
bys_model <- naiveBayes(class~. ,data = train[,-c(1,2)])
##贝叶斯预测
bys_pred <- predict(bys_model,newdata = test[,-c(1,2)])
##提取分类概率
bys_prob <- predict(bys_model,newdata = test[,-c(1,2)],type = "raw")
##生成四格表
bys_perf<-table(test$class,bys_pred,dnn=c("Actual","Predicted"))
bys_perf
##查看
performance(bys_perf)
##ROC
bys_roc <- roc(test$class,bys_prob[,2])
plot(bys_roc,print.auc=T)
plot(smooth(bys_roc),print.auc=T)

##逻辑回归
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
seed <- 37
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
log_model <- glm(class~.,data = train[,-c(1,2)],family = "binomial")
log_prob <- predict(object =log_model,newdata=test[,-c(1,2)],type = "response")
log_pred<-ifelse(log_prob>=0.5,"synergistic","antagonistic")
log_pred <- factor(log_pred)
##ROC
log_roc <- roc(test$class,log_prob)
plot(log_roc,print.auc=T)
plot(smooth(log_roc),print.auc=T)

##ANN
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
seed <- 150
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)

##模型
ann_model <- nnet(class~.,data = train[,-c(1,2)],size=10)
#测试集
ann_pred <- predict(ann_model,newdata = test[,-c(1,2)],type="class")
ann_prob <- predict(ann_model,newdata = test[,-c(1,2)])
###ROC曲线
ann_roc <- roc(test$class,ann_prob[,1])
plot(ann_roc, print.auc=TRUE)
plot(smooth(ann_roc),print.auc=T)

##投票整合部分
setwd("github\\NetPro\\Data\\")
##设置种子
seed <- 431
set.seed(seed)
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
seed <- 431
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##数据
set <- rbind(pos,neg)

##SVM
seed <- 3
set.seed(seed)
svm_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
svm_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(svm_train)))),]
set.seed(seed)
svm_model <- svm(class~.,data = svm_train[,-c(1,2)],kernel="radial",
                 cost=1,gamma=1/ncol(svm_train[,-c(1,2)]),probability = TRUE)
svm_pred <- predict(svm_model,newdata = svm_test[,-c(1,2)])

##Logistic分类器
seed <- 3
set.seed(seed)
log_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
log_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(log_train)))),]
set.seed(seed)
log_model <- glm(class~.,data = log_train[,-c(1,2)],family = "binomial")
log_prob <- predict(object =log_model,newdata=log_test[,-c(1,2)],type = "response")
log_pred<-ifelse(log_prob>=0.5,"synergistic","antagonistic")
log_pred <- factor(log_pred)

##ANN分类器
seed <- 6
set.seed(seed)
ann_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
ann_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(ann_train)))),]
set.seed(seed)
ann_model <- nnet(class~.,data = ann_train[,-c(1,2)],size=10)
ann_pred <- predict(ann_model,newdata = ann_test[,-c(1,2)],type="class")
ann_pred <- factor(ann_pred)
names(ann_pred) <- rownames(ann_test)

##NB分类器
seed <- 6
set.seed(seed)
bys_train <- set[unique(sample(nrow(set),nrow(set),replace = T)),]
bys_test <- set[which(as.numeric(rownames(set))%in%setdiff(as.numeric(rownames(set)),as.numeric(rownames(bys_train)))),]
set.seed(seed)
bys_model <- naiveBayes(class~. ,data = bys_train[,-c(1,2)])
bys_pred <- predict(bys_model,newdata = bys_test[,-c(1,2)])
names(bys_pred) <- rownames(bys_test)

##交集信息
test_ins <- c(rownames(svm_test),rownames(bys_test),
              rownames(log_test),rownames(ann_test))
test_ins <- as.data.frame(table(test_ins))
test_ins$test_ins <- as.character(test_ins$test_ins)
test_ins <- test_ins[test_ins$Freq>=3,]

##添加类别
for (i in 1:nrow(test_ins)) {
  if(length(svm_pred[names(svm_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$SVM[i] <- svm_pred[names(svm_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$SVM[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(bys_pred[names(bys_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$Bys[i] <- bys_pred[names(bys_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$Bys[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(log_pred[names(log_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$Log[i] <- log_pred[names(log_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$Log[i] <- NA
  
}

for (i in 1:nrow(test_ins)) {
  if(length(ann_pred[names(ann_pred)%in%test_ins$test_ins[i]])!=0)
  {test_ins$ANN[i] <- ann_pred[names(ann_pred)%in%test_ins$test_ins[i]]
  }else
    test_ins$ANN[i] <- NA
  
}
##去除频率信息
test_ins <- test_ins[,-2]
##投票判断
for (i in 1:nrow(test_ins)) {
  if(length(which(test_ins[i,-1]%in%2))>=3){
    test_ins$pred_class[i] <- "synergistic"
  }else
    test_ins$pred_class[i] <- "antagonistic"
}
##添加正式类别
for (i in 1:nrow(test_ins)) {
  test_ins$actual_class[i] <- as.character(set$class[which(rownames(set)%in%test_ins$test_ins[i])])
}
test_ins$pred_class <- factor(test_ins$pred_class)
test_ins$actual_class <- factor(test_ins$actual_class)
##整合ROC曲线
int_roc <- roc(as.numeric(test_ins$actual_class),
               as.numeric(test_ins$pred_class))
plot(int_roc,print.auc=T)

##将ROC曲线进行整合
##ggplot2绘图
roc <- ggroc(list("RF (AUC=0.879)"=rf_roc,
                  "ANN + LR + SVM + NB (AUC=0.779)"=int_roc,
                  "SVM (AUC=0.764)"=svm_roc ,
                  "NB (AUC=0.775)"=bys_roc,
                  "LR (AUC=0.762)"=log_roc,
                  "ANN (AUC=0.749)"=ann_roc),legacy.axes = T,
             size=1)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 16),#设置y轴标题的字体属性
        axis.title.x=element_text(size = 16),#设置x轴标题的字体属性
        panel.border = element_rect(),
        axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        #legend.key.height  = unit(30, "pt"),
        #legend.key.width = unit(30, "pt"),
        legend.background    =element_rect(color="black"),
        legend.position = c(.7,.18),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,1,1),'pt'))+
  ylab("Sensitivity")+xlab("1 - Specificity")
roc <- roc+theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))
roc

#########################################
#3.Use RF model to view the AUC for different measures
##读入药物组合金标准
setwd("github\\NetPro\\Data\\")
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)

##保留六个测度
seed <- 431
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,-c(1,2)],mtry=2,
                         ntree=1000,
                         important=TRUE,proximity=TRUE)

##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,-c(1,2)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,-c(1,2)],"prob")
rf_roc <- roc(test$class,rf_prob[,2])
plot(rf_roc,print.auc=T)

##6个测度整合模型
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
##删除PCA和RWR
std <- std[,-c(9,10)]
##保留六个测度
seed <- 1
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,-c(1,2)],mtry=2,
                         ntree=1000,
                         important=TRUE,proximity=TRUE)

##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,-c(1,2)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,-c(1,2)],"prob")
rf_roc6 <- roc(test$class,rf_prob[,2])
plot(rf_roc6,print.auc=T)

##结构
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
std <- std[,c(1,2,3,11)]
seed=89
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

##设置种子训练结构模型
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,c(3,11)],
                         ntree=1000,
                         important=TRUE,proximity=TRUE)
##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,c(3,11)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,c(3,11)],"prob")
##ROC曲线
stru_roc <- roc(test$class,rf_prob[,2])
plot(stru_roc,print.auc=T)

##ATC
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
std <- std[,c(1,2,4,11)]
seed=62
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

##设置种子训练结构模型
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,c(3,4)],
                         ntree=1000,
                         important=TRUE,proximity=TRUE)
##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,c(3,4)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,c(3,4)],"prob")
##生成四格表
rf_perf<-table(test$class,rf_pred,dnn=c("Actual","Predicted"))
##ROC曲线
atc_roc <- roc(test$class,rf_prob[,2])
plot(atc_roc,print.auc=T)

##Dis
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
std <- std[,c(1,2,5,11)]
seed=36
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

##设置种子训练结构模型
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,c(3,4)],
                         ntree=1000,
                         important=TRUE,proximity=TRUE)
##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,c(3,4)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,c(3,4)],"prob")
##生成四格表
rf_perf<-table(test$class,rf_pred,dnn=c("Actual","Predicted"))
##ROC曲线
dis_roc <- roc(test$class,rf_prob[,2])
plot(dis_roc,print.auc=T)

##Jaccard
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
std <- std[,c(1,2,6,11)]
seed=22
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

##设置种子训练结构模型
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,c(3,4)],
                         ntree=1000,
                         important=TRUE,proximity=TRUE)
##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,c(3,4)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,c(3,4)],"prob")
##生成四格表
rf_perf<-table(test$class,rf_pred,dnn=c("Actual","Predicted"))
##ROC曲线
tar_roc <- roc(test$class,rf_prob[,2])
plot(tar_roc,print.auc=T)

##GO
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
std <- std[,c(1,2,7,11)]
seed=51
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

##设置种子训练结构模型
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,c(3,4)],
                         ntree=1000,
                         important=TRUE,proximity=TRUE)
##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,c(3,4)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,c(3,4)],"prob")
##ROC曲线
go_roc <- roc(test$class,rf_prob[,2])
plot(go_roc,print.auc=T)

##KEGG
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
std <- std[,c(1,2,8,11)]
seed=92
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

##设置种子训练结构模型
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,c(3,4)],
                         ntree=1000,
                         important=TRUE,proximity=TRUE)
##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,c(3,4)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,c(3,4)],"prob")
##ROC曲线
kegg_roc <- roc(test$class,rf_prob[,2])
plot(kegg_roc,print.auc=T)

##PCA
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
std <- std[,c(1,2,9,11)]
seed=27
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

##设置种子训练结构模型
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,c(3,4)],
                         ntree=1000,
                         important=TRUE,proximity=TRUE)
##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,c(3,4)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,c(3,4)],"prob")
##ROC曲线
pca_roc <- roc(test$class,rf_prob[,2])
plot(pca_roc,print.auc=T)

##RWR
std <- as.data.frame(fread("Gold standard for drug combination.csv"))
std$class <- factor(std$class)
std <- std[,c(1,2,10,11)]
seed=18
set.seed(seed)
#区分阴性和阳性
pos <- std[std$class=="synergistic",]
neg <- std[std$class=="antagonistic",]
##设置种子
set.seed(seed)
##抽取1：1
neg <- neg[sample(nrow(neg),nrow(pos)),]
##设置种子
set.seed(seed)
##抽取train和test
train <- rbind(pos[sample(nrow(pos),0.7*nrow(pos)),],neg[sample(nrow(neg),0.7*nrow(neg)),])
std <- rbind(pos,neg)
test <- std[which(as.numeric(rownames(std))%in%setdiff(as.numeric(rownames(std)),as.numeric(rownames(train)))),]

##设置种子训练结构模型
set.seed(seed)
rf_model <- randomForest(class ~ ., data=train[,c(3,4)],
                         ntree=1000,
                         important=TRUE,proximity=TRUE)
##随机森林预测
rf_pred <- predict(rf_model,newdata=test[,c(3,4)])#"prob"显示概率
##概率
rf_prob <- predict(rf_model,newdata=test[,c(3,4)],"prob")
##ROC曲线
rwr_roc <- roc(test$class,rf_prob[,2])
plot(rwr_roc,print.auc=T)

##ROC曲线图
p_roc <- ggroc(list("NetPro-based (AUC=0.88)"=rf_roc,
                    "Six SM-based (AUC=0.81)"=rf_roc6,
                    "CSS-based (AUC=0.69)"=stru_roc,
                    "ATCS-based (AUC=0.77)"=atc_roc,
                    "TNDS-based (AUC=0.74)"=dis_roc,
                    "TSOS-based (AUC=0.73)"=tar_roc,
                    "GOS-based (AUC=0.67)"=go_roc,
                    "DIPS-based (AUC=0.69)"=kegg_roc,
                    "DDS-based (AUC=0.61)"=pca_roc,
                    "RWRS-based (AUC=0.61)"=rwr_roc),legacy.axes = T,
               size=1)+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,colour = "black"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(size=12,colour = "black"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(size = 16),#设置y轴标题的字体属性
        axis.title.x=element_text(size = 16),#设置x轴标题的字体属性
        panel.border = element_rect(size = 1),
        axis.line = element_line(colour = "black"), #去除默认填充的灰色，并将x=0轴和y=0轴加粗显示(size=1)
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        #legend.key.height  = unit(0.6, "cm"),
        #legend.key.width = unit(1, "cm"),
        legend.background    =element_rect(color="black"),
        legend.position = c(.73,.32),
        panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),
        plot.margin=unit(c(1,1,1,1),'pt'),
        plot.title = element_text(hjust = 0.5,size = 16))+
  ylab("Sensitivity")+xlab("1 - Specificity")+ggtitle("ROC")
p_roc
p_roc <- p_roc+theme(plot.margin = unit(c(0.7,1,0.7,0.8), "cm"))
p_roc