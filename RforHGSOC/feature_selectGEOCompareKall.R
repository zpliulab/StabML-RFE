## clear
rm(list = ls())

## package
library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用


## External input data
setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataGEOTogather")
# gene = as.matrix(read.csv("DataPythonK\\ dif_feature_xgbrfe .csv", header=TRUE, sep = ',')[1:14,])

# Internal input data
fren = read.csv("fren.csv", header=TRUE, sep = ',')
gene <- as.matrix(fren[which(fren[,2] >= 2),1])

# feature extraction ------------------------------------------------------

## name function
library(stringr)
myname <- function(x){
  # x <- 1
  name <- dir[x]
  name1 <- str_split_fixed(name, "./", 2);
  name2 <- str_split_fixed(name1[2], "/", 2)
  name3 <- str_split_fixed(name2[2], "[_]", 2);
  name6 <- name3[1]
  return(name6)
}

## clcle 
myfile = list.files("DataGEO")    
dir = paste("./DataGEO/", myfile, sep="")   
n = length(dir) 


for (i in 1:n) {
  # i <- 1
  Data1 <- read.table(file = dir[i], header = T, check.names = FALSE)
  colnames(gene) <- c('gene')
  # Data1 <- t(Data)
  Data2 <- cbind(rownames(Data1), Data1)
  colnames(Data2) <- c('gene', colnames(Data1))
  # View(Data2[,1:10])
  
  genedata <- merge(gene, Data2, by = "gene")#[,-2]
  genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(Data1[1,],genedata1)
  rownames(genedata2) <- c('Lable', rownames(genedata1))
  
  path <- paste("./DataCompare/GEO/",paste(str_c(myname(i),"_bio"),".txt"))
  write.table(genedata2, path, quote=F,sep="\t")
  
  
  data1 = read.table("GSE6942827651_outcome.txt", header = T, check.names = FALSE)
  dim(data1)
  gene <- as.matrix(genedata[,1])
  
  colnames(gene) <- c('gene')
  # Data1 <- t(Data)
  data2 <- cbind(rownames(data1), data1)
  colnames(data2) <- c('gene', colnames(data1))
  # View(data2[,1:10])
  
  genedata <- merge(gene, data2, by = "gene")#[,-2]
  genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
  genedata2 <- rbind(data1[1,],genedata1)
  rownames(genedata2) <- c('Lable', rownames(genedata1))
  path <- paste("./DataCompare/TCGA/",paste(str_c("TCGA_",myname(i),"_bio"),".txt"))
  write.table(genedata2, path, quote=F, sep="\t")
}



# classo performance  -----------------------------------------------------

## clear
rm(list = ls())

## package
library(pROC)
library(e1071)


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

## input data
setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataGEOTogather\\DataCompare")


myfile = list.files("GEO")    
dir = paste("GEO/", myfile, sep="")   
n = length(dir) 


for (i in 1:n) {
  # i <- 5
  myfile = list.files("GEO")    
  dir = paste("GEO/", myfile, sep="") 
  Data  = read.table(file = dir[i], header = T, check.names = FALSE)
  
  myfile = list.files("TCGA")    
  dir = paste("TCGA/", myfile, sep="") 
  Data2 = read.table(file = dir[i], header = T, check.names = FALSE)
  
  x.train <- data.frame(t(Data2)[,-1])
  y.train <- t(Data2)[,1]
  x.test <- data.frame(t(Data)[,-1])
  y.test <- t(Data)[,1]
  
  # model <- svm(x.train, y.train, kernel = 'linear', scale = FALSE)   # linear  radial
  # summary(model)
  
  model <- glm(y.train~., data = x.train, family = binomial)
  # glmfit <- glm(y.train~., data = x.train, family = binomial, control = list(maxit = 100))
  summary(model)
  
  ## prediction
  p_test <- predict(model, x.test, type = "response")
  p_test = as.matrix(p_test)
  A_test <- data.frame(p_test, y.test)
  names(A_test)<- c("p", "outcome")
  
  # p <- A_test[,1]
  # p_glm <- cbind(log(p/(1-p)), A_test[,2])
  # colnames(p_glm) <- c('y.test', 'Lable')
  # write.table(p_glm,"D:\\E\\博士\\R_程序\\HGSOC\\DataPRS\\pre40595cvlog.txt",quote=F,sep="\t")

  
  path <- paste("./DataAUC/",paste(myname(i),".csv"))
  write.csv(A_test, path, row.names = F)
  index <- classindex(A_test, t(x.test), 0.500)
  picauc <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
  print(picauc$auc)
  print(index)
}





