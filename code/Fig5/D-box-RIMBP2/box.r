library(ggplot2)
# setwd('E:\\mine\\work\\pancreatic\\version3\\code\\Fig5\\D-box-RIMBP2')
data=read.table('RIMBP2.txt',sep='\t',header=TRUE,row.names=1)

p<-ggplot(data,aes(ATAC1,log(ATAC+1,2)))+labs( y= 'ATAC log2(openness+1)')+geom_boxplot( fill="white") +geom_jitter(width = 0.2)+
	#scale_fill_brewer(palette="Dark2") +
	theme_classic()+theme(legend.position='bottom',axis.text.x = element_text())
pdf("RIMBP2.pdf",width=2,height=3)
grid.draw(p)
dev.off()

p<-ggplot(data,aes(RNA1,log(RIMBP2+1,2)))+labs(y= 'RNA log2(FPKM+1)')+geom_boxplot( fill="white") +geom_jitter(width = 0.2)+
	#scale_fill_brewer(palette="Dark2") +
	theme_classic()+theme(legend.position='bottom',axis.text.x = element_text())
pdf("RIMBP2-RNA.pdf",width=2,height=3)
grid.draw(p)
dev.off()


