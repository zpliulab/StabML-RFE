## clear
rm(list = ls())

## package
library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用


## input data
setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataTcgaGtex")

## internal validation
# gene = as.matrix(read.csv("DataPythonK\\ dif_feature_gbmrfe .csv", header=TRUE, sep = ',')[1:16,])

## external validation
fren = read.csv("fren4.csv", header=TRUE, sep = ',')
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
  write.table(genedata2, path, ,quote=F,sep="\t")
  
  
  data1 = read.table("matrix_DEoutcome.txt", header = T, check.names = FALSE)
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

## input data
setwd("D:\\E\\博士\\R_程序\\HGSOC\\Datatcgagtex\\DataCompare")


myfile = list.files("GEO")    
dir = paste("GEO/", myfile, sep="")   
n = length(dir) 


for (i in 1:n) {
  # i <- 2
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
  
  # model <- svm(x.train, y.train, kernel = 'linear', scale = FALSE)   # linear   radial
  # summary(model)
  # p_test <- predict(model, x.test, type = "response")
  
  # glm.fit <- glm(y.train~., data = x.train, family = binomial, control = list(maxit = 100))
  glm.fit <- glm(y.train~., data = x.train, family = binomial)
  summary(glm.fit)
  # glm.fit$coefficients
  p_test <- predict(glm.fit, x.test, type = "response")
  
  
  p_test = as.matrix(p_test)
  A_test <- data.frame(p_test, y.test)
  names(A_test)<- c("p", "outcome")
  # write.csv(A_test,"A_test_54388.csv",row.names = F)
  
  
  # p <- A_test[,1]
  # p_glm <- cbind(log(p/(1-p)), A_test[,2])
  # colnames(p_glm) <- c('y.test', 'Lable')
  # setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataPRS")
  # write.table(p_glm,"pre59491_73685_4cvlog.txt",quote=F,sep="\t")
  
  # pdf(file = "ROC_GSE27651_52.pdf",width = 5,height = 5)
  show_25 <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
  print(i)
  print(show_25$auc)
}
