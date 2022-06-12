

##2022.6.6 SPTB
rm(list = ls())

## 
library(EFS)


setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataTcgaGtex\\DataPythonKauc")
dataAB <- as.matrix(read.csv("dif_feature_ABrfe .csv", header = T))
dataNB <- as.matrix(read.csv("dif_feature_NBrfe.csv", header = T))
# dataNNET <- as.matrix(read.csv("dif_feature_NNETrfe .csv", header = T))
dataSVM <- as.matrix(read.csv("dif_feature_SVMrfe .csv", header = T))



# uniongene <- as.matrix(union(dataAB, union(dataNB, union(dataNNET, dataSVM))))
uniongene <- as.matrix(union(dataAB, union(dataNB, dataSVM)))
uniongene <- cbind(uniongene, rep(0, dim(uniongene)[1]))
colnames(uniongene) <- c("union", "rank")


rank <- rep(27:1,1)/27
dataABrank <- cbind(dataAB, rank)
colnames(dataABrank) <- c("AB-RFE", "rank")
# dataNNETrank <- cbind(dataNNET, rank)
# colnames(dataNNETrank) <- c("NNET-RFE", "rank")
dataSVMrank <- cbind(dataSVM, rank)
colnames(dataSVMrank) <- c("SVM-RFE", "rank")
dataNBrank <- cbind(dataNB, rank)
colnames(dataNBrank) <- c("NB-RFE", "rank")



weight <- function(uniongene,dataABrank){
  for (i in 1:dim(uniongene)[1]) {
    for (j in 1:dim(dataABrank)[1]) {
      if(uniongene[i,1]==dataABrank[j,1]){
        uniongene[i,2]=dataABrank[j,2]
      }
    }
  }
  colnames(uniongene) <- colnames(dataABrank)
  return(uniongene)
}

a <- uniongene
dataga <- as.matrix(rbind(weight(a, dataABrank)[,2], weight(a, dataNBrank)[,2], weight(a, dataSVMrank)[,2]))
colnames(dataga) <- a[,1]
dataga[1,2]
dataga <- apply(dataga, 2, as.numeric)
dataga[1,2]
rownames(dataga) <- c("AB-RFE","NB-RFE","SVM-RFE")
dataga[1,2]


number = 3
table = dataga/number


### part3 
# barplot_fs("test", efs, order = TRUE)

name <- "SPTB"
efs_table <- table
order = TRUE


if(order == TRUE){
  efs_table <- efs_table[, order(colSums(efs_table))]
}

paranr = length(efs_table[1,])  
if(paranr>100){
  b =  colSums(efs_table)
  #b= a[order(a)]
  
  pdf(paste(name,'.pdf', sep=""),
      width= 12,
      height= 12)
  barplot(b,
          ylim=c(0,1),
          main= 'Ensemble Feature Selection',
          xlab = "Features",
          ylab = "Importance values",
          axisnames = FALSE
  )
  dev.off()
}
if(paranr<35){h=10}else {h=(paranr/5)}

names=colnames(efs_table)
cols = c('goldenrod1','navy',
         'royalblue','indianred3',
         'darkolivegreen1','darkgreen',
         'darkolivegreen3','chartreuse4')
# pdf(paste(name,'.pdf', sep=""),
# width= 12,
# height= h)
par(mar=c(5, 4, 4, 10), xpd=TRUE)
barplot= barplot(efs_table,
                 xlim=c(0,0.7),
                 # xlim=c(0,1),
                 main= 'StabML-RFE feature selection for SPTB',
                 horiz=T,
                 # las=2,  # 横坐标竖着 90 度
                 las=1,
                 names.arg=abbreviate(names),
                 # names.arg=names,
                 col=cols)
legend("topright", inset=c(-0.2,0), legend=row.names(efs_table),col=cols, lty=1, lwd=12 )
text(colSums(efs_table)+0.035,barplot,
     format(round(colSums(efs_table), 2),T))
segments(1, 0, 1, 1.25*paranr, lty = 3, col = "gray40")
# dev.off()


