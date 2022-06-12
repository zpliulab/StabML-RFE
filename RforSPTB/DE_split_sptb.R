## clear
rm(list=ls())


library(dplyr)
library(tidyverse)



## makedir
setwd("D:\\E\\²©Ê¿\\R_³ÌÐò\\HGSOC\\DataTcgaGtex")

# data split --------------------------------------------------------------

datascale <- as.matrix(read.table("matrix_DEsptb.txt",header=T,sep='\t', check.names = F))
label <- rbind(as.matrix(rep(c("1"), 98)), as.matrix(rep(c("0"), 228)))
colnames(label) <- c("outcome")
datalab <- rbind(t(label), datascale)

# datalab <- rbind(label, datascale)
# write.table(datalab, file = "matrix_DEoutcome.txt",quote = F, sep = "\t")
datalabt <- data.frame(t(datalab))


set.seed(2022)

library(caret)
library(dplyr)
trainingsamples <- datalabt$outcome %>% createDataPartition(p = 0.7, list = FALSE)
datanum <- datascale
traindata  <- datanum[, trainingsamples]
testdata <- datanum[, -trainingsamples]
View(traindata[,1:10])
View(testdata[,1:10])
traindata[100,3]

# write.table(traindata, file = "matrix_DEtrain.txt",quote = F, sep = "\t")
# write.table(testdata, file = "matrix_DEtest.txt",quote = F, sep = "\t")



label <- rbind(as.matrix(rep(c("1"), 98)), as.matrix(rep(c("0"), 228)))
colnames(label) <- c("outcome")
trainLab <- as.numeric(label[trainingsamples])
testLab <- as.numeric(label[-trainingsamples])
# View(trainLab)
sum(trainLab)

# write.table(trainLab, file = "trainLab.txt",quote = F, sep = "\t")
# write.table(testLab, file = "testLab.txt",quote = F, sep = "\t")
