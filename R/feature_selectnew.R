## clear
rm(list = ls())

## package
library(dplyr)        
library(tidyr)
library(tidyverse)    

## input data
setwd("D:\\E\\博士\\R_程序\\HGSOC")
Data1 = read.table("HGSOC的参考数据\\GSE120196\\GSE120196_scale.txt", header = T, check.names = FALSE)
gene = read.csv("datanew\\bio.csv", header=TRUE, sep = ',')
dim(Data1)

colnames(gene) <- c('gene')
Data2 <- cbind(rownames(Data1), Data1)
colnames(Data2) <- c('gene', colnames(Data1))

genedata <- merge(gene, Data2, by = "gene")#[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(Data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
write.table(genedata2,"datanew\\GSE120196_bio.txt",quote=F,sep="\t")

# 在原数据集提取特征 ---------------------------------------------------------------
data1 = read.table("data\\GSE69428_scale.txt", header = T, check.names = FALSE)
dim(data1)
gene <- as.matrix(genedata[,1])

colnames(gene) <- c('gene')
data2 <- cbind(rownames(data1), data1)
colnames(data2) <- c('gene', colnames(data1))
genedata <- merge(gene, data2, by = "gene")#[,-2]
genedata1 <- genedata %>% tibble::column_to_rownames(colnames(.)[1])
genedata2 <- rbind(data1[1,],genedata1)
rownames(genedata2) <- c('Lable', rownames(genedata1))
write.table(genedata2,"datanew\\GSE69428_120196_bio.txt",quote=F,sep="\t")
