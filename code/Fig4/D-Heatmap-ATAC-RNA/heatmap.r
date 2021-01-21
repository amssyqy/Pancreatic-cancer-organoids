# setwd('E:\\mine\\work\\pancreatic\\version3\\code\\Fig4\\D-Heatmap-ATAC-RNA')
library(pheatmap)
EXP_c=read.table('final-network-max-exp.txt',sep='\t',header=FALSE)
ATAC_c=read.table('final-network-max-enh.txt',sep='\t',header=FALSE)
ATAC_c=as.matrix(ATAC_c)
ATAC_zscore=ATAC_c
for (i in 1:dim(ATAC_c)[1])
{
  mean0=mean(ATAC_c[i,])
  sigma=sd(ATAC_c[i,])
  ATAC_zscore[i,]=(ATAC_c[i,]-mean0)/sigma
}
EXP_c=as.matrix(EXP_c)
EXP_zscore=EXP_c
for (i in 1:dim(EXP_c)[1])
{
  mean0=mean(EXP_c[i,])
  sigma=sd(EXP_c[i,])
  EXP_zscore[i,]=(EXP_c[i,]-mean0)/sigma
}
EXP_zscore[EXP_zscore>2]=2
EXP_zscore[EXP_zscore<=-2]=-2
ATAC_zscore[ATAC_zscore>2]=2
ATAC_zscore[ATAC_zscore<=-2]=-2
color1=colorRampPalette(c("#3363aa","#1384d5","#23b2ae","#c5bc5e","#faf513"))(100)
library(pheatmap)
p=pheatmap(EXP_zscore,color=color1,show_rownames=FALSE,treeheight_row=0,treeheight_col=0)
RNA_1=EXP_zscore[,p$tree_col$order]
ATAC_1=ATAC_zscore[,p$tree_col$order]


RNA1_promoter=RNA_1[1:365,]
ATAC1_promoter=ATAC_1[1:365,]
RNA1_pos=RNA_1[366:(365+3185),]
ATAC1_pos=ATAC_1[366:(365+3185),]
RNA1_neg=RNA_1[(365+3185):(365+3185+2072),]
ATAC1_neg=ATAC_1[(365+3185):(365+3185+2072),]
p=pheatmap(RNA1_promoter,cluster_col=FALSE,color=color1,show_rownames=FALSE,treeheight_row=0,treeheight_col=0)
RNA1_promoter=RNA1_promoter[p$tree_row$order,]
ATAC1_promoter=ATAC1_promoter[p$tree_row$order,]
p=pheatmap(RNA1_pos,cluster_col=FALSE,color=color1,show_rownames=FALSE,treeheight_row=0,treeheight_col=0)
RNA1_pos=RNA1_pos[p$tree_row$order,]
ATAC1_pos=ATAC1_pos[p$tree_row$order,]
p=pheatmap(RNA1_neg,cluster_col=FALSE,color=color1,show_rownames=FALSE,treeheight_row=0,treeheight_col=0)
RNA1_neg=RNA1_neg[p$tree_row$order,]
ATAC1_neg=ATAC1_neg[p$tree_row$order,]

RNA=rbind(RNA1_promoter,RNA1_pos,RNA1_neg)
ATAC=rbind(ATAC1_promoter,ATAC1_pos,ATAC1_neg)


color0=colorRampPalette(c("#3363aa","#18b3fd","#b4d3e0","#fab672","#b41f1f"))(100)
library(ggplot2)


p1=pheatmap(ATAC,show_rownames=FALSE,cluster_col=FALSE,cluster_row=FALSE,color=color0,show_colnames=FALSE,treeheight_row=0,treeheight_col=0,gaps_row=c(365,365+3185))
pdf('ATAC-enhancer.pdf',width=4,height=10)
p1
dev.off()


p=pheatmap(RNA,color=color1,cluster_col=FALSE,cluster_row=FALSE,show_rownames=FALSE,show_colnames=FALSE,treeheight_row=0,gaps_row=c(365,365+3185),treeheight_col=0)
pdf('EXP-promoter.pdf',width=4,height=10)
p
dev.off()
