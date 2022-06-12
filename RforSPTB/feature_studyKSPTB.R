rm(list = ls())

library(dplyr)       # ％>％ 管道函数的调用，传参
library(tidyr)
library(tidyverse)   # tibble 的调用



setwd('D:\\E\\博士\\R_程序\\HGSOC\\DataTcgaGtex')

# 提取 feature gene 表达值 -----------------------------------------------------------------
Data1 = read.table("matrix_DEoutcome.txt", header = T, check.names = FALSE)
# View(Data1[,1:10])
fren = read.csv("fren4.csv", header=TRUE, sep = ',')

## 2022.6.8
# fren = read.csv("fren5.csv", header=TRUE, sep = ',')
# View(fren)

## frenquence
gene <- as.matrix(fren[which(fren[,2] >= 2),1])
colnames(gene) <- c('gene')
# write.csv(gene, file = "bio.csv", row.names = F)

# Data1 <- t(Data)
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))
# View(Data2[,1:10])

genedata <- merge(gene, Data2, by = "gene")#[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
dim(genedata2)    #  19 20
write.table(genedata2,"datainter\\GSE59491_bio.txt",quote=F,sep="\t") 

