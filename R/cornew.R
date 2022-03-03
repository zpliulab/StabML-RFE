
## clear 
rm(list=ls())

## load data
setwd("D:\\E\\博士\\R_程序\\HGSOC")

x0 = read.table("Datapython20\\GSE69428_bio.txt", header = T, sep='\t', fill=TRUE, strip.white = T, check.names=F)
rownames(x0)
x1 <- data.frame(t(x0))
colnames(x1) <- rownames(x0)
x2 <- x1[,-1]

# heatmap -----------------------------------------------------------------
library(pheatmap)

selected <- t(x2)
Label = read.table("data\\phe_zh.txt",header = TRUE, sep = "\t")
Label <- factor(Label[,"class"])
Label <- data.frame(Label)
rownames(Label) = colnames(selected)

p <- pheatmap(selected[1:18,],annotation_col = Label, 
              color = colorRampPalette(c("blue", "white","red"))(100),
              fontsize_row = 8,scale = "row", cutree_cols = 2, border_color = NA)
p
# pdf(file = "datapython20\\PicHeatmap.pdf",width = 8,height = 8)
p
# dev.off()



# DE gene bar -------------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(dplyr) 
library(data.table)
library(ggsignif)

x = read.table("Datanew\\GSE69428_bio_forDEGs.txt", header = T, sep='\t', fill=TRUE, strip.white = T)
pred_log <- data.frame(t(x))
colnames(pred_log) <- rownames(x)
dim(pred_log)

# Compare two independent groups ------------------------------------------
k = 1
p <- ggboxplot(pred_log, x = "Lable", y = gene[k],
               color = "Lable", palette = "jco",
               add = "jitter")
p + stat_compare_means(method = "t.test")


# add p value -------------------------------------------------------------
annotation=c("1.9e-05", "3.4e-05", "1.6e-08", "4.4e-06", "1.4e-06",
             "1.5e-05", "2.8e-07", "8.1e-10", "3.6e-07", "1.2e-06",
             "8.1e-07", "4.8e-05", "1.4e-05", "4.7e-05", "9.4e-08",
             "6.7e-08", "1.8e-08", "3.8e-06")



# 变成长格式的数据： ---------------------------------------------------------------

b <- reshape2::melt(pred_log, id.vars = c("Lable"))
b$Label <- as.factor(b$Lable)


## 做每个箱线图的均值点及均值连线，需获得每个组每个属性的均值，
## 并定义每个组每个属性的X坐标为固定值 
# group1的mean:
c<-copy(b)
setDF(c)
c1<-tapply(c[c$Label==1,"value"],c[c$Label==1,"variable"],mean)
c2<-tapply(c[c$Label==0,"value"],c[c$Label==0,"variable"],mean)
c3<-rbind(data.frame(variable=names(c1),value=c1,Label=1),data.frame(variable=names(c2),value=c2,Label=0))

c3$Label<-as.factor(c3$Label)
## average
c3$variable2<-NA

gene <- c3$variable 
lab1 <- c()
lab2 <- c()

## lab1,lab2 position
lab1[1] = 1.2
lab2[1] = 0.8
c3[c3$Label==1&c3$variable==gene[1],"variable2"]<-lab1[1]
c3[c3$Label==0&c3$variable==gene[1],"variable2"]<-lab2[1]

## 固定显著性的位置
y <- c()
y[1] = max(max(pred_log[1:10,gene[1]]),max(pred_log[11:20,gene[1]]) ) +0.5

for (l in 2:18) {
  y[l] = max(max(pred_log[1:10,gene[l]]),max(pred_log[11:20,gene[l]]) ) +0.5
}

data=data.frame(x=lab2,
                xend=lab1,
                y=y,
                annotation=annotation)


for (i in 2:18) {
  # i = 2
  lab1[i] <- lab1[i-1] + 1
  lab2[i] <- lab2[i-1] + 1
  c3[c3$Label==1&c3$variable==gene[i],"variable2"]<-lab1[i]
  c3[c3$Label==0&c3$variable==gene[i],"variable2"]<-lab2[i]
}

p1<-ggplot(b)+
  geom_boxplot(aes(x=variable,y=value,fill=Label),width=0.6,
               position = position_dodge(0.8),outlier.size = 0,outlier.color = "white")+
  scale_fill_manual(values = c("red", "blue"),breaks=c("1","0"),labels=c("Normal","HGSOC"))+
  geom_point(data=c3,aes(x=variable2,y=value,color=Label),shape=15,size=1)+
  geom_line(data=c3,aes(x=variable2,y=value,color=Label),size=1,linetype = "dotted")+
  # geom_smooth(data=c3,method = 'loess',formula = 'y ~ x',
  #             aes(x=variable2,y=value,color=Label),size=1,linetype = "dashed")+
  # stat_summary(fun = mean, geom = "errorbar",
  #              aes(x=variable,y=value,ymax = ..y.., ymin = ..y..,color=Label),
  #              width = .75, linetype = "dashed") +
  xlab("")+
  ylab("")+
  scale_y_continuous(limits = c(1, 14),breaks=seq(1, 14, 2)) +
  geom_signif(stat="identity",
              data=data.frame(x=lab2,
                              xend=lab1,
                              y=y,
                              annotation=annotation),
              aes(x=x,xend=xend, y=y, yend=y, 
                  annotation=annotation, hjust=0, vjust=0.5, angle=90)) +
  xlab("Biomarker") +
  ylab("Expression value") +
  theme_bw()+
  theme(
    legend.position = "right",
    # legend.background=element_blank(),
    # legend.key = element_blank(),
    # legend.margin=margin(0s,0,0,0,"mm"),
    axis.text.x = element_text(size=rel(1.1),angle=25, colour = "black", vjust = 0.7), # face="bold",
    axis.line.x = element_line(size = 0.5, colour = "black"),
    axis.line.y = element_line(size = 0.5, colour = "black"),
    legend.text=element_text(size=rel(1.1)),
    legend.title=element_blank()
    # panel.border = element_blank(),
    # panel.grid = element_blank()
  ) # +
  #guides(color=FALSE)
p1

# pdf(file = "datanew\\PicDEGs11.pdf",width = 11.5,height = 4)
p1
# dev.off()


# corralation -------------------------------------------------------------
library(Hmisc) 
dataposet <- x2
res2 <- rcorr(as.matrix(dataposet))    # pearson(默认) spearman
View(res2$r[1:10,1:10])
View(res2$P[1:10,1:10])

l <- 18
library(PerformanceAnalytics)
chart.Correlation(dataposet[,1:l], histogram=TRUE, pch=19, fontsize_row = 10)
p2 <- chart.Correlation(dataposet[,1:l], histogram=TRUE, pch=19)
# pdf(file = "datanew\\Picchart.pdf",width = 10,height = 10)
p2
# dev.off()

library(corrplot)
data0 <- res2$r[1:l,1:l]

cols<-c('#3E5CC5','#65B48E','#E6EB00','#E64E00')
pal <- colorRampPalette(cols)
image(x=1:l,y=1,z=as.matrix(1:l),col=pal(l))

p3 <- corrplot(data0, type="lower", 
         col = pal(10), #c("purple", "green","white"),
         tl.cex = 0.6,  #指定文本标签的大小 
         tl.col = "black",  #指定文本标签的颜色 
         tl.srt = 45
) 
# pdf(file = "datanew\\Piccor.pdf")
# p3
# dev.off()


# cor net -----------------------------------------------------------------
cor_pearson <- cor(x2, method = 'pearson')
class(cor_pearson)
# fix(cor_pearson)
cor_pearson[abs(cor_pearson) <= 0.8] <- 0
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
# write.csv(ed,"datanew\\cor_select08_83.csv",quote = F,row.names = F)

# 提取gene,得到logFC ------------------------------------------------------------------
##############  手动将. 换成-  #####################
ed = read.csv("datanew\\cor_select08_83.csv", header=TRUE, sep = ',')
gene1 <- ed[,1]
gene2 <- ed[,2]
cor_0.9_gene <- as.matrix(union(gene1,gene2) ) # 26
rownames(cor_0.9_gene) <- cor_0.9_gene[,1]

# dif文件
dif_FC = read.csv("data\\dif_FC.csv", header=TRUE, sep = ',')
rownames(dif_FC) <- dif_FC[,1]
dif_FC_gene <- dif_FC[rownames(cor_0.9_gene),]

# 下调
gene_neg <- dif_FC_gene[which(dif_FC_gene$logFC < 0),1]
write.csv(gene_neg,"datanew\\gene_neg.csv",quote = F,row.names = F)
# 上调
gene_pos <- dif_FC_gene[which(dif_FC_gene$logFC > 0),1]
write.csv(gene_pos,"datanew\\gene_pos.csv",quote = F,row.names = F)


cor_select2 <- list()
cor_select1 <- as.matrix(cor_select)
for (i in 1 : dim(ed)[1]){
  if (ed[i,1] == cor_select1[i,1] & ed[i,2] == cor_select1[i,2]){
    cor_select2 <- rbind(cor_select2, cor_select1[i,])
  }
}


g <- p1
myc <- clusters(g, mode="strong")
mycolor <- c('green', 'yellow', 'red', 'skyblue')
V(g)$color <- mycolor[myc$membership + 1]
plot(g, layout=layout.fruchterman.reingold,  
     vertex.size=4,  
     vertex.label = V(p1)$name,  
     vertex.label.cex=0.7,  
     vertex.label.dist=0.4, 
     vertex.label.color = "black"  
)



