## clear
rm(list = ls())

## package
library(pROC)
library(ggplot2)
library(ROCR)  

## data
setwd("D:\\E\\²©Ê¿\\R_³ÌÐò\\HGSOC\\Datanew")

roc1 <- read.csv("A_test_26712.csv", header=TRUE, sep = ',')
roc2 <- read.csv("A_test_40595.csv", header=TRUE, sep = ',')
roc3 <- read.csv("A_test_120196.csv", header=TRUE, sep = ',')
roc4 <- read.csv("A_test_27651.csv", header=TRUE, sep = ',')
roc5 <- read.csv("A_test_54388.csv", header=TRUE, sep = ',')

myroc <- function(matrix){
  pred <- prediction(matrix[,1], matrix[,2])  
  perf <- performance(pred,"tpr","fpr") 
  x <- unlist(perf@x.values)   
  y <- unlist(perf@y.values)
  plotdata <- data.frame(x,y) 
  names(plotdata) <- c("x", "y")
  return(plotdata)
}
plotdata <- myroc(roc1)
ggplot(plotdata) + 
  geom_path(aes(x = x, y = y, colour = x), size=1) + 
  labs(x = "False positive rate", y = "Ture positive rate") +     
  scale_colour_gradient(name = 'False positive rate', low = 'blue', high = 'red') +
  # theme(plot.title = element_text(face = 'bold',size=15))
  theme_bw() + coord_equal()


Roc1 <- roc(roc1[,2],roc1[,1])
g <- ggroc(Roc1)
g
g + theme_minimal() +  
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="red", 
               linetype=6) 
gl <- ggroc(Roc1, legacy.axes = TRUE)
gl

Roc2 <- roc(roc2[,2],roc2[,1])
Roc3 <- roc(roc3[,2],roc3[,1])
Roc4 <- roc(roc4[,2],roc4[,1])
Roc5 <- roc(roc5[,2],roc5[,1])

g2 <- ggroc(list(GSE26712=Roc1, 
                 GSE40595=Roc2, 
                 GSE120196=Roc3, 
                 GSE27651=Roc4,
                 GSE54388=Roc5),
            legacy.axes = TRUE)
g2

g2 + annotate(geom = "segment", 
              x = 0, y = 0, xend =1, yend = 1, 
              colour = "gray", size = 0.5) +
  scale_fill_discrete(labels=c("GSE26712", "GSE40595", "GSE120196", "GSE27651", "GSE54388")) +
  theme_gray() + coord_equal() +
  theme(legend.position = c(0.70,0.25), #legend.position = 'inside',
        legend.text = element_text(color = 'black',size = 12),
        axis.text = element_text(color = 'black',size = 15),
        axis.text.x = element_text(angle = 0),
        axis.title = element_text(color = 'black',size = 15, family = 'Arial', face = 'plain'),
        axis.ticks = element_line(color = 'black'))
		