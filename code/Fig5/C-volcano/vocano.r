library(ggplot2)
library(ggthemes)
# library(Cairo)
library(ggrepel)
# setwd('E:\\mine\\work\\pancreatic\\version3\\code\\Fig5\\C-volcano')
data=read.table("forplot.txt", sep='\t',head=T)

### 火山图
data$threshold = 'Unchanged'
data$threshold[(data$FC>1.5)]='Mut up'
data$threshold[(data$FC<1/1.5)]='Mut down'
genename=read.table('gene-20201129.txt',header=TRUE)

ll=as.character(genename$gene)
ll[ll=='0']=''
p=ggplot(data=data, aes(x=log(FC,2), y=opennes, colour=threshold, fill=threshold)) + 
 scale_color_manual(values=c('Mut up'="purple", 'Unchanged'="grey",'Mut down'="orange"))+
 geom_point(size=0.8) +geom_text_repel(aes(label = ll), size = 2)+
 xlim(c(-4, 4)) +
 ylim(c(0, 50)) +
 theme_bw() +
 geom_vline(xintercept=c(-log(1.5,2),log(1.5,2)),lty=4,lwd=0.6)+
 theme(legend.position="right",
 panel.grid=element_blank(),
 legend.title = element_blank(),
 legend.text= element_text(face="bold", color="black", size=8),
 plot.title = element_text(hjust = 0.5),
 axis.text.x = element_text(face="bold", color="black", size=8),
 axis.text.y = element_text(face="bold", color="black", size=8),
 axis.title.x = element_text(face="bold", color="black", size=8),
 axis.title.y = element_text(face="bold",color="black", size=8))+
 labs(x="log2 (Fold Change)",y="Average Openness)")
 library(grid)
pdf('vocano.pdf',width=4,height=3)
grid.draw(p)
dev.off()



