## clear
rm(list = ls())

## package
library(pROC)
library(e1071)

## input data
setwd("D:\\E\\博士\\R_程序\\HGSOC\\Datanew")

# GSE73685 7个gene (17个样本)-----------------------------------------------------------------
Data  = read.table("GSE54388_bio.txt", header = T, check.names = FALSE)
Data2 = read.table("GSE69428_54388_bio.txt", header = T, check.names = FALSE)

# 独立数据集 -------------------------------------------------------------------
x.train <- data.frame(t(Data2)[,-1])
y.train <- t(Data2)[,1]
x.test <- data.frame(t(Data)[,-1])
y.test <- t(Data)[,1]

# SVM ---------------------------------------------------------------------
## 使用svm函数训练支持向量机
set.seed(666) 
tuned <- tune.svm(x.train,y.train, gamma = 10^(-6:-1), cost = 10^(1:2)) # tune
summary (tuned) # to select best gamma and cost

model <- svm(x.train, y.train, kernel = "radial", cost = tuned$best.parameters$cost, gamma=tuned$best.parameters$gamma,  scale = FALSE)
summary(model)

## 利用svm() 函数建立的模型进行预测
p_test <- predict(model, x.test, type = "response")
p_test = as.matrix(p_test)
A_test <- data.frame(p_test, y.test)
names(A_test)<- c("p", "outcome")
write.csv(A_test,"A_test_54388.csv",row.names = F)
classindex(A_test, t(Data), 0.836)

# pdf(file = "ROC_GSE27651_52.pdf",width = 5,height = 5)
show_25 <- plot.roc(A_test$outcome, A_test$p, print.auc=T, main="pAUC")
show_25$auc
# s1 <- smooth(show_25,method="binormal")
# plot(s1)

## 分类指标的函数
classindex <- function(A_test, data, thro){
  predict = ifelse(A_test[,1] >= thro, 1, 0)
  predict_value = predict
  true_value = A_test[,2]
  error = predict_value-true_value
  
  data <- t(Data)
  # 计算模型准确性（accuracy）、精确度（Precision），召回率（Recall）-- 敏感性（sensitivity）-- 真阳性率（TPR）和 F测度（F-measure）和混淆矩阵
  accuracy = (nrow(data)-sum(abs(error)))/nrow(data)
  precision = sum(true_value & predict_value)/sum(predict_value)  
  recall = sum(predict_value & true_value)/sum(true_value)        
  F_measure = 2*precision*recall/(precision+recall)   
  specificity = ((nrow(data)-sum(abs(error)))-sum(predict_value & true_value))/(nrow(data)-sum(true_value))
  table = table(predict_value, true_value) 
  result <- list(accuracy, precision, recall, specificity, F_measure,table)
  names(result) <- c("accuracy", "precision", "recall", "specificity", "F_measure", "table")
  return(result)
}

