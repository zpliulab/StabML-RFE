
rm(list = ls())
setwd("D:\\E\\博士\\R_程序\\ABvsCD\\DataNew")


# bar plot --------------------------------------------------------------------
setwd("D:\\E\\博士\\R_程序\\HGSOC\\Datanew")
data0 <- read.csv("result_ROC_appro.csv", header=F, sep = ',')[,1:6]
data1 <- t(data0)[-1,]
colnames(data1) <- c("Index", "GSE26712", "GSE40595", "GSE120196", "GSE27651", "GSE54388")
data2 <- data.frame(data1)
data <- data2


library(reshape2)     
library(ggplot2)

mydata <- melt(data, id.vars="Index", 
               variable.name="Dataset",
               value.name="Value")
ggplot(mydata, aes(Index, Value, fill = Dataset)) + 
  geom_bar(stat="identity", position="dodge")


## change color
library(ggplot2)
library(ggthemes)

sub <- factor(mydata$Index,levels = c('Acc', 'Pre', 'Sn', 'Sp', 'F_measure'))

ggplot(mydata, aes(Dataset, Value, fill = Index)) + 
  # geom_bar(stat="identity", position="dodge") +
  geom_bar(aes(fill = sub), stat="identity", position="dodge", width=.5) +
  theme_hc() + 
  ## theme_bw(),theme_minimal(), theme_classic(), theme_gray(), theme_excel(), 
  ## theme_economist(), theme_fivethirtyeight(), theme_tufte(), theme_geocs(), 
  ## theme_calc(), theme_hc()
  scale_fill_wsj("colors6", "") +     
  #scale_fill_wsj(palette = "colors6", ...)       
  #palette    character The color palette to use: . "rgby", "red_green", "black_green", "dem_rep", "colors6"
  guides(fill=guide_legend(title=NULL)) +
  # ggtitle("The number of all metabolite and differential metabolite") + 
  theme(axis.title = element_text(color = 'black',size = 14),
        axis.text = element_text(color = 'black',size = 12), # element_blank()
        # axis.text = element_text(angle = 45, color = 'black',size = 12, family = 'Arial', face = 'plain'), # element_blank()
        legend.position = 'top', # 'none'  c(0.9,1.2)
        legend.title = element_blank(),
        legend.text = element_text(color = 'black',size = 14),
        axis.text.x = element_text(angle = 0))
