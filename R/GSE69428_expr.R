
## clear 
rm(list=ls())

## packages
library(stringr)
library(BiocGenerics)
library(parallel)
library(Biobase)
library(dplyr)       
library(tidyr)
library(tidyverse)   
library(fdrtool)      
library(data.table)   
library(fdrtool)      
library(data.table)   

## load data
setwd("D:\\E\\博士\\R_程序\\HGSOC\\Data")
exprSet <- read.table("GSE69428_series_matrix.txt",header=T,sep='\t',fill=TRUE,strip.white = T)
exprSet$ID_REF <- as.character(exprSet$ID_REF)


## load annoation file
anno <- read.table("GPL570.txt",header=T,sep='\t',fill=TRUE,strip.white = T,quote = "")   
anno2 <- anno %>%      
  select(ID,Gene.ID) %>%             
  filter(Gene.ID != '') 
colnames(anno2) <- c('ID_REF','EntrezID')    
anno2$ID_REF <- as.character(anno2$ID_REF)   


## 将基因表达数据与芯片注释文件的探针名进行对应
exprset2 <- exprSet %>%                      
  inner_join(anno2,by='ID_REF') %>%         
  select(ID_REF,EntrezID, everything())   
                 

## 整理芯片注释文件，把其中一个探针对应多个基因的拆分开
exprset3 <- exprset2
a <- tibble(exprset3[,1:2])
test1 <- apply(a,1, function(x){
  str_split(x[2],'///',simplify=T)     
} )


test2 <- apply(a, 1, function(x){         
  paste(x[1],str_split(x[2],'///', simplify=T), sep = "---")
})


unlist(test2)                             
x <- tibble(unlist(test2))                
colnames(x) <- "lala"                      
x2 <- separate(x,lala,c("id","entrezID"),sep = '---')     
x2[1:10,1]
expset3[1:10,1]
x3 <- merge(x2,exprset3,by.x = "id",by.y="ID_REF",all=FALSE)   
x4<-x3[,-c(1,3)]                          


zz <- as.matrix(apply(as.matrix(x4[,1]),1,function(x) as.numeric(x)))
View(zz)


XX <- x4[,-1]
colnames(XX)[1:3]
XX1 <- cbind(zz,XX)
colnames(XX1) <- c("entrezID",colnames(XX))


## 用基因id对整理好的芯片注释文件进行基因名的更新
homo<-read.table("homo.txt",header=T,sep='\t')
homo[1:6,1]
x5 <- merge(homo,XX1,by.x="GeneID",by.y = "entrezID",all=FALSE) 


## 探针名匹配基因名，取出多个探针对应一个基因的数据计算IQR，保留IQR最大的探针数据
expset4 <- x5 %>%
  dplyr::select(-GeneID) %>%              
  mutate(rowIQR =apply(.[,-1],1,IQR)) %>%  
  arrange(desc(rowIQR)) %>%                
  distinct(Symbol,.keep_all = T) %>%       
  dplyr::select(-rowIQR) %>%                  
  tibble::column_to_rownames(colnames(.)[1])  
View(expset4[1:10,1:10])
dim(expset4)   # 17689    29
# write.table(expset4,"GSE69428_expr.txt",quote=F,sep="\t")  


## load label
lable = read.csv("GSE69428_all.csv", header = T, sep=',')
lable_1 <- lable[which(lable$Sample.type == "tissue"),]
expset5 <- expset4[,which(lable$Sample.type == "tissue")]
dim(expset5)          # 17689    20  


lable2 = read.csv("GSE69428_20.csv", header = T, sep=',')
dim(lable2)    # 20   2
data = rbind(as.matrix(t(lable2[,2])), as.matrix(expset5))
rownames(data) <- c('Lable', rownames(expset5))
View(data[1:10,])
dim(data)    # 17690    20
# write.table(data,"GSE69428_outcome.txt",quote=F,sep="\t")   



## scale
data = read.table("GSE69428_outcome.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names = F)

my_scale <- function(x){
  x1 <- cbind(t(x[1,]), scale(t(x[-1,])))
  x2 <- t(x1)
  return(x2)
}

data1 <- my_scale(data)
dim(data1)    # 17690    20
# write.table(data1,"GSE69428_scale.txt",quote=F,sep="\t")  

