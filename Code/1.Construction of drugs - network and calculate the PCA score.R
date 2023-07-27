#######github项目
##1.Construct drug-drug network

##载入包
library(data.table)
library(psych)

##读入不同药物测度的矩阵信息
setwd("github\\NetPro\\Data\\")
##读入整合文件
com <- as.data.frame(fread("feature_matrix.txt"))
colnames(com) <- c("Drug1","Drug2","DDS","ATC","Structure","GO","KEGG","Jaccard","Distance")
com <- com[,-3]
##删除NA，保证每个测度都有信息
com <- na.omit(com)

##只保留单一测度的信息
dt <- com[,3:8]
##计算相关性矩阵
sigma <- cor(dt)
##计算特征值和特征向量
e <- eigen(sigma)
##特征值
e$values
##特征向量
e$vectors
#提取结果中的特征值，即各主成分的方差；
val <- e$values
##查看大于1的方差个数
which(val>1 )
#换算成标准差(Standard deviation);
std_sev <- sqrt(val)
#计算方差贡献率和累积贡献率
pov <- val/sum(val)
cp <- cumsum(pov)

##画碎石图,确定主成分个数
fa.parallel(dt,fa="pc",show.legend = F)
abline(h=1)
##碎石图方法2
par(mar=c(6,6,2,2))
plot(e$values,type="b",
     cex=1,
     cex.lab=2,
     cex.axis=2,
     lty=1,
     lwd=1,
     col="blue",
     pch=4,
     xlab = "Component Number",
     ylab="eigen values of principal components")
abline(h=1,col="red",lty=3)

##提取主成分
fit <- principal(dt,nfactors = 2,rotate = "varimax",#方差最大旋转
                 scores = T)
fit

##绘制主成分分析荷载矩阵
fa.diagram(fit,digits = 2,rsize = 0.8)

##读入fit中测度在不同主成分的系数和不同主成分比例系数
fit
##读入手动写出fit中的系数
w <- as.data.frame(fread("Measure coefficient.csv"))
rownames(w) <- w$V1
w <- w[,-1]
##主成分权重
w_PCA <- as.data.frame(fread("Principal component coefficient.csv"))

##使用权重计算得分
for (i in 1:nrow(com)) {
  com$RC1[i] <- sum(com[i,3:8]*w$RC1)
  com$RC2[i] <- sum(com[i,3:8]*w$RC2)
  com$PCA[i] <- w_PCA$RC1*com$RC1[i]+w_PCA$RC2*com$RC2[i]
  print(i)
}

##写出文件
write.table(com,"Drug measure with PCA score.csv",sep = ",",col.names = T,row.names = F)
