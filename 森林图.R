library(grid)
library(forestploter)

setwd("D:/双细胞孟德尔随机化/133_双细胞孟德尔随机化/double_cell/1_bc91")
mydata <- read.table("森林图.txt",header = T,sep = "\t")
mydata$pvalue=ifelse(mydata$pvalue<0.001, "<0.001", sprintf("%.4f", mydata$pvalue))
mydata$` ` <- paste(rep(" ", 30), collapse = " ")
mydata$'OR(95%CI)'<-ifelse(is.na(mydata$or),"",
                           sprintf('%.4f(%.4f to %.4f)',
                                   mydata$or,mydata$or_lci,mydata$or_uci))
mydata[is.na(mydata)] <- " "

mydata$"P value"=mydata$pvalue
tm <- forest_theme(base_size = 8,
                   ci_pch = 20,
                   ci_col = "#4575b4",
                   ci_lty = 1,
                   ci_lwd = 2.3,
                   ci_Theight = 0.2, 
                   refline_lwd = 1.5,
                   refline_lty = "dashed",
                   refline_col = "red",
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   footnote_cex = 1.1,
                   footnote_fontface = "italic",
                   footnote_col = "blue")
forest(mydata[,c(1:3, 10:12)],
       est = mydata$or,
       lower = mydata$or_lci,
       upper = mydata$or_uci,
       sizes = 0.6,
       ci_column = 4,
       ref_line = 1,
       xlim = c(0,2),
       ticks_at = c(0,1,2),
       xlab="OR",
       theme = tm)