data0 = read.table('2000gene-exp-f.txt',header=TRUE,sep='\t',row.names=1)
library(ggplot2)
library(ggrepel)
library(ggbiplot)
library(RColorBrewer)
library(grid)
##color for tables
cols <- c("PDAC" = "grey", "NEN" = "#e64b35", "ACC" = "#3c5488", "NPO" = "#91d1c2","IPMN"='#ffe100')
mydata1 = t(as.matrix(data0))
mydata.pca <- prcomp(log(mydata1+1))
scores <- mydata.pca$x
scores=data.frame(scores)
rownames(scores)=colnames(data0)
##labels
scores$Subtype = c("PDAC",'IPMN',rep("PDAC",33),"NEN","NEN","ACC",rep("PDAC",3),'NEN',rep("PDAC",2),rep("NPO",3),'NEN')
N1=as.character(round(summary(mydata.pca)$importance[3,1],2))
N2=as.character(round(summary(mydata.pca)$importance[3,2]-summary(mydata.pca)$importance[3,1],2))
N3=as.character(round(summary(mydata.pca)$importance[3,3]-summary(mydata.pca)$importance[3,2],2))
xl=paste("PC1 ","(",N1,")",sep='')
yl=paste("PC2 ","(",N2,")",sep='')
zl=paste("PC3 ","(",N3,")",sep='')
p <- ggplot(scores, aes(PC1, PC2,colour = Subtype)) +
    geom_point(size=3) +
    xlab(xl) + ylab(yl)  + scale_colour_manual(values = cols) +
    theme(panel.background = element_rect(fill = "white", colour = "black"),axis.text=element_text(size=7),axis.title=element_text(size=7,face="bold"))
pdf(file='PCA-log_point.pdf',width= 4,height=2.5)
grid.draw(p)
dev.off()


