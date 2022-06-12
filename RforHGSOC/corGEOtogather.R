
## clear
rm(list = ls())

setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataGEOtogather")

## load data
x0 = read.table("DataInter\\GSE6942827651_bio.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names=F)
rownames(x0)
x1 <- data.frame(t(x0))
colnames(x1) <- rownames(x0)
x2 <- x1[,-1]


label <- cbind(rownames(x2),x1[,1])
label[which(label[,2] == "1"),2] <- c("Normal")
label[which(label[,2] == "0"),2] <- c("HGSOC")
colnames(label) <- c("sample","class")
# write.table(label, file = "phe_HGSOCGEOtogather.txt",quote=F,sep="\t",row.names = F)

# heatmap -----------------------------------------------------------------
library(pheatmap)

selected <- t(x2)
Label = read.table("phe_HGSOCGEOtogather.txt",header = TRUE, sep = "\t")
Label <- factor(Label[,"class"])
Label <- data.frame(Label)
rownames(Label) = colnames(selected)

p <- pheatmap(selected,annotation_col = Label, 
              color = colorRampPalette(c("blue", "white","red"))(100),
              fontsize_row = 8,scale = "row", cutree_cols = 2, border_color = NA)
p



# Box DEGs ----------------------------------------------------------------

# https://www.zhihu.com/question/44891604

library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(cowplot)


Exp_plot <- x2
Exp_plot$sam <- Label$Label
Exp_plot$sam <- factor(Exp_plot$sam,levels=c("HGSOC","Normal"))


col <-c("#F0AD4E","#D9534F")  # "#5CB85C","#337AB7","#F0AD4E","#D9534F"

gene <- read.csv("bio.csv",header=T)
gene[2,1]


plist2 <- list()#创建一个空列表，用来存储循环的产出
for (i in 1:(dim(gene)[1])){
  # i <- dim(gene)[1]
  bar_tmp<-Exp_plot[,c(gene[i,1],"sam")]#循环提取每个基因表达信息
  colnames(bar_tmp)<-c("Expression","sam")#统一命名
  my_comparisons1 <- list(c("HGSOC", "Normal")) #设置比较组
  # my_comparisons2 <- list(c("Asymptomatic", "Severe"))#设置比较组
  # my_comparisons3 <- list(c("Asymptomatic", "Critical"))#设置比较组
  # my_comparisons4 <- list(c("Mild", "Severe"))#设置比较组
  # my_comparisons5 <- list(c("Mild", "Critical"))#设置比较组
  # my_comparisons6 <- list(c("Severe", "Critical"))#设置比较组
  pb1<-ggboxplot(bar_tmp,#ggboxplot画箱线图
                 x="sam",#x轴为组别
                 y="Expression",#y轴为表达量
                 color="sam",#用样本分组填充
                 fill=NULL,
                 add = "jitter",#添加散点
                 bxp.errorbar.width = 0.6,
                 width = 0.4,
                 size=0.01,
                 font.label = list(size=30), 
                 palette = col)+theme(panel.background =element_blank())
  pb1<-pb1+theme(axis.line=element_line(colour="black"))+theme(axis.title.x = element_blank())#坐标轴修饰
  pb1<-pb1+theme(axis.title.y = element_blank())+theme(axis.text.x = element_text(size = 15,angle = 0,vjust = 1,hjust = 0.5))#横坐标文字设置
  pb1<-pb1+theme(axis.text.y = element_text(size = 15))+ggtitle(gene[i,1])+theme(plot.title = element_text(hjust = 0.5,size=15,face="bold"))#标题设置
  pb1<-pb1+theme(legend.position = "NA")#（因为有组图，横坐标分组了，所以不需要设置legend）
  pb1<-pb1+stat_compare_means(method="t.test",hide.ns = F,
                              comparisons =c(my_comparisons1),
                              # comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3,my_comparisons4,my_comparisons5,my_comparisons6),
                              label="p.signif")  # p.format  p.signif
  #显著性检验用t检验，添加不同比较组。详情可以查看stat_compare_means函数帮助信息
  plist2[[i]]<-pb1 #将画好的图储存于plist2列表，并不断赋值循环直到结束
}

plot_grid(plist2[[1]],plist2[[2]],plist2[[3]],
          plist2[[4]],plist2[[5]],plist2[[6]],
          plist2[[7]],plist2[[8]],plist2[[9]],
          plist2[[10]],plist2[[11]],plist2[[12]],
          plist2[[13]],plist2[[14]],plist2[[15]],
          plist2[[6]],plist2[[7]],plist2[[8]],
          plist2[[9]],plist2[[10]])

## delete the later five pictures


# corralation -------------------------------------------------------------

# 计算相关系数及显著性
library(Hmisc) 
dataposet <- x2
res2 <- rcorr(as.matrix(dataposet))    # pearson(默认) spearman


l <- 15
library(PerformanceAnalytics)
chart.Correlation(dataposet[,1:l], histogram=TRUE, pch=19, fontsize_row = 10)
p2 <- chart.Correlation(dataposet[,1:l], histogram=TRUE, pch=19)
# pdf(file = "Picchart.pdf",width = 15,height = 15)
p2
# dev.off()

library(corrplot)
data0 <- res2$r[1:l,1:l]

cols<-c('#3E5CC5','#65B48E','#E6EB00','#E64E00')
pal <- colorRampPalette(cols)
image(x=1:l,y=1,z=as.matrix(1:l),col=pal(l))

p3 <- corrplot(data0, type="lower", 
         col = pal(10), #c("purple", "green","white"),
         tl.cex = 0.8,  #指定文本标签的大小 
         tl.col = "black",  #指定文本标签的颜色 
         tl.srt = 45
) 
# pdf(file = "Piccor.pdf")
p3
# dev.off()


# cor net -----------------------------------------------------------------
cor_pearson <- cor(x2, method = 'pearson')
class(cor_pearson)
# fix(cor_pearson)
#提取“cor_pearson”中相关系数 >0.5 或 <-0.5 的结果输出为 csv 样式
cor_pearson[abs(cor_pearson) <= 0.60] <- 0
# write.csv(cor_pearson, file='cor_pearson.csv', quote = FALSE)

library(modelr)
library(dbplyr)
library(tidyverse)

correlation <- cor_pearson
cor <- correlation %>% as.tibble(.) %>%
  add_column(Var1=rownames(correlation)) %>%
  gather(rownames(correlation),key='Var2',value='cor')
cor_select <- cor[cor$cor != 0 & cor$cor != 1,]

library(igraph)
PP <- graph_from_data_frame(cor_select,directed = F)
p1 <- simplify(PP,remove.loops = T,remove.multiple = T)  # 最终的数对
ed <- as_edgelist(p1, names = TRUE)
union(ed[,1], ed[,2])
# write.csv(ed,"cor_select06_42.csv",quote = F,row.names = F)


# 提取gene,得到logFC ------------------------------------------------------------------
##############  手动将. 换成-  #####################

ed = read.csv("cor_select06_42.csv", header=TRUE, sep = ',')
gene1 <- ed[,1]
gene2 <- ed[,2]
cor_0.9_gene <- as.matrix(union(gene1,gene2) ) # 26
rownames(cor_0.9_gene) <- cor_0.9_gene[,1]

## dif文件
dif_FC = read.csv("dif_FC.csv", header=TRUE, sep = ',')
rownames(dif_FC) <- dif_FC[,1]
dif_FC_gene <- dif_FC[rownames(cor_0.9_gene),]

# ## 下调
gene_neg <- dif_FC_gene[which(dif_FC_gene$logFC < 0),1]
# View(gene_neg)
# write.csv(gene_neg,"gene_neg.csv",quote = F,row.names = F)
# ## 上调
gene_pos <- dif_FC_gene[which(dif_FC_gene$logFC > 0),1]
# View(gene_pos)
# write.csv(gene_pos,"gene_pos.csv",quote = F,row.names = F)



