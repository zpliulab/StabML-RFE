
## clear
rm(list = ls())

# load(".RData")    # 载入工作区

library(e1071)
library(limma)
library(ggplot2)
library(reshape2)
library(gmodels)
library(pROC)

# 读入数据 --------------------------------------------------------------------
setwd("D:\\E\\博士\\R_程序\\HGSOC\\Datatcgagtex")
eset <- as.matrix(read.table("matrix_DEtrain.txt",header = TRUE,sep = "\t"))
# View(eset[,1:10])

## name function
library(stringr)
myname <- function(x){
  # x <- 1
  name <- dir[x]
  name1 <- str_split_fixed(name, "./", 2);
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[_]", 2);
  name4 <- str_split_fixed(name3[2], "[.]", 2);
  name6 <- str_c("dif_feature_",name4[1])
  return(name6)
}


## 分类指标的函数
classindex <- function(A_test, Data, thro){
  predict = ifelse(A_test[,1] >= thro, 1, 0)
  predict_value = predict
  true_value = A_test[,2]
  error = predict_value-true_value
  
  data <- t(Data)
  # 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
  accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
  precision = sum(true_value & predict_value)/sum(predict_value)  #真实值预测值全为1 / 预测值全为1 --- 提取出的正确信息条数/提取出的信息条数
  recall = sum(predict_value & true_value)/sum(true_value)        #真实值预测值全为1 / 真实值全为1 --- 提取出的正确信息条数 /样本中的信息条数
  # P和R指标有时候会出现的矛盾的情况，这样就需要综合考虑他们，最常见的方法就是F-Measure（又称为F-Score）
  F_measure = 2*precision*recall/(precision+recall)    #F-Measure是Precision和Recall加权调和平均，是一个综合评价指标
  specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))
  table = table(predict_value, true_value) 
  result <- list(accuracy, precision, recall, specificity, F_measure,table)
  names(result) <- c("accuracy", "precision", "recall", "specificity", "F_measure", "table")
  return(result)
}

## clcle 
myfile = list.files("DataPython")    
dir = paste("./DataPython/", myfile, sep="")   
n = length(dir) 


for (i in 1:n) {
  # i <- 1
  rfe_svm <- read.table(file = dir[i], stringsAsFactors = FALSE)
  colnames(rfe_svm) <- c("rank","chipID")
  ## rate
  alpha <- 0.045
  k <- floor((dim(eset)[1])*alpha)    # 向下取整
  rfe_svm_dif <- head(rfe_svm, n=k)    # n = 40L
  eset_svm <- eset[rfe_svm_dif$chipID,]
  path <- paste("./DataPythonK/",paste(myname(i),".csv"))
  write.csv(rfe_svm_dif[,2], path, row.names = F)
  
  
  
  ## predict data
  x <- as.data.frame(t(eset_svm))
  y <- read.table("trainLab.txt",header = TRUE,sep = "\t")[,1]
  test <- read.table("matrix_DEtest.txt",header = TRUE,sep = "\t")
  xtest <- as.data.frame(t(test[rfe_svm_dif$chipID,]))
  ytest <- read.table("testLab.txt",header = TRUE,sep = "\t")[,1]
  
  
  ## train
  set.seed(666) 
  # tuned <- tune.svm(x, y, gamma = 10^(-6:-1), cost = 10^(1:2)) # tune
  # summary (tuned) # to select best gamma and cost
  # model <- svm(x, y, kernel = 'linear')
  # # model <- svm(x, y1, kernel = "radial", cost = tuned$best.parameters$cost, gamma=tuned$best.parameters$gamma,  scale = FALSE)
  # summary(model)
  
  model <- glm(y~., data = x, family = binomial)
  # model <- glm(y~., data = x, family = binomial, control = list(maxit = 100))
  summary(model)


  ## prediction
  p_test <- predict(model, xtest, type = "response")
  p_test = as.matrix(p_test)
  A_test <- data.frame(p_test, ytest)
  names(A_test)<- c("p", "outcome")
  path <- paste("./DataAUC/",paste(myname(i),".csv"))
  write.csv(A_test, path, row.names = F)
  index <- classindex(A_test, t(xtest), 0.500)
  picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
  print(picauc$auc)
  print(index)
}






