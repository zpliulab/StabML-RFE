library(ggplot2)
library(ggpubr)
library(Rcpp)

pred_log = read.table("D:\\E\\博士\\R_程序\\HGSOC\\DataPRS\\pre59491_73685_4cvlog.txt", header = T, check.names = FALSE)

dim(pred_log)

term0 <- which(pred_log$Lable == 0)
pred_log[term0, 2] <- "Term"
term1 <- which(pred_log$Lable == 1)
pred_log[term1, 2] <- "SPTB"

colnames(pred_log) <- c("preterm_risk_score","label")

# View(pred_log)
head(pred_log)
class(pred_log)
###########################################   画图   ##############################################

pred_log$lable <- as.factor(pred_log$label)  # dose 为列名，即表达值

e <- ggplot(pred_log, aes(x = label, y = preterm_risk_score))  # x 有三类，就有3个箱子

# Basic box plot
e + geom_boxplot() 


e + geom_boxplot(aes(fill = label), position=position_dodge(0.8)) +
  geom_jitter(aes(color = label, shape = label),
              position=position_dodge(0.8))

 test --------------------------------------------------------------------

ee <- ggboxplot(pred_log, x = "label", y = "preterm_risk_score", color = "label",
                palette = "npg", add = "jitter") 

my_comparisons <- list( c("SPTB", "Term") )
# my_comparisons <- list( c("Treatment", "Control") )
ee + stat_compare_means(comparisons = my_comparisons,method = "wilcox.test")  # t.test


ee + stat_compare_means(comparisons = my_comparisons) + 
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "blue")


setwd("D:\\E\\博士\\R_程序\\HGSOC\\DataPRS")
# jpeg(file = "Box_59491_73685_4cvlog.jpg")
ee + stat_compare_means(comparisons = my_comparisons) + 
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "blue") 
# + stat_compare_means(method = "t.test", label.x = 1.38, label.y = 2.4)
# dev.off()

# over --------------------------------------------------------------------


# Box plot with mean points
e + geom_boxplot() +
  stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 4, color = "blue")
# Change box plot colors by groups
e + geom_boxplot(aes(fill = label))

# mean points and color
e + geom_boxplot(aes(fill = label))+
  stat_summary(fun.y = mean, geom = "point",
               shape = 18, size = 4, color = "blue")

## 保存

# jpeg(file = "Box_59491_73685_4cvlog.jpg")
# pdf(file = "Box_59491_73685_4cvlog.pdf")
e + geom_boxplot(aes(fill = lable))+
  stat_summary(fun = mean, geom = "point",
               shape = 18, size = 4, color = "blue")
# dev.off()

##########################################################################################
