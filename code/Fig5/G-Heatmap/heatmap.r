# setwd('E:\\mine\\work\\pancreatic\\version3\\code\\Fig5\\G-Heatmap')
D=read.table('drug-20200709.txt',sep='\t',header=TRUE, row.names=1)
ATAC=read.table('merge_region_head_quantiler-rerank.txt',sep='\t',header=TRUE, row.names=1)
D=as.matrix(D)
ATAC1=as.matrix(log(ATAC+1,2))
N=dim(D)[1]
M=dim(ATAC)[1]
# A=cor(t(ATAC1),t(D))
A=cor(t(ATAC1),t(D),method='pearson')

choosed_ATAC=c()
choosed_drug=c()
corr0=c()
choosed_ATAC=c()
choosed_drug=c()

for (i in 1:M)
{
for (j in 1:N)
{
if (A[i,j]>0.5)
{
choosed_ATAC=c(choosed_ATAC,i)
choosed_drug=c(choosed_drug,j)
}
}
}
drug1=D[choosed_drug,]
ATAC_1=ATAC1[choosed_ATAC,]
corr0=c()
for (i in 1:length(choosed_ATAC))
{
corr0=c(corr0,cor.test(ATAC_1[i,],drug1[i,])$p.value)
}

choosed=(corr0<0.01)
ATAC_c=ATAC_1[choosed,]
EXP_c=drug1[choosed,]


ATAC_zscore=ATAC_c
for (i in 1:dim(ATAC_c)[1])
{
  mean0=mean(ATAC_c[i,])
  sigma=sd(ATAC_c[i,])
  ATAC_zscore[i,]=(ATAC_c[i,]-mean0)/sigma
}
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


choosed_ATAC=c()
choosed_drug=c()

for (i in 1:M)
{
for (j in 1:N)
{
if (A[i,j]<=-0.5)
{
choosed_ATAC=c(choosed_ATAC,i)
choosed_drug=c(choosed_drug,j)
}
}
}
drug1=D[choosed_drug,]
ATAC_1=ATAC1[choosed_ATAC,]
corr0=c()
for (i in 1:length(choosed_ATAC))
{
corr0=c(corr0,cor.test(ATAC_1[i,],drug1[i,])$p.value)
}

choosed=(corr0<0.01)
ATAC_c=ATAC_1[choosed,]
EXP_c=drug1[choosed,]


ATAC_zscore1=ATAC_c
for (i in 1:dim(ATAC_c)[1])
{
  mean0=mean(ATAC_c[i,])
  sigma=sd(ATAC_c[i,])
  ATAC_zscore1[i,]=(ATAC_c[i,]-mean0)/sigma
}
EXP_zscore1=EXP_c
for (i in 1:dim(EXP_c)[1])
{
  mean0=mean(EXP_c[i,])
  sigma=sd(EXP_c[i,])
  EXP_zscore1[i,]=(EXP_c[i,]-mean0)/sigma
}
EXP_zscore1[EXP_zscore1>2]=2
EXP_zscore1[EXP_zscore1<=-2]=-2
ATAC_zscore1[ATAC_zscore1>2]=2
ATAC_zscore1[ATAC_zscore1<=-2]=-2

EXP_z=rbind(EXP_zscore,EXP_zscore1)
ATAC_z=rbind(ATAC_zscore,ATAC_zscore1)
color0=colorRampPalette(c("#3363aa","#18b3fd","#b4d3e0","#fab672","#b41f1f"))(100)

library(pheatmap)
library(ggplot2)

p=pheatmap(EXP_z)

EXP_1=EXP_z[,p$tree_col$order]
ATAC_1=ATAC_z[,p$tree_col$order]
EXP1=EXP_1[1:11084,]
ATAC1=ATAC_1[1:11084,]
EXP2=EXP_1[11085:(11084+4313),]
ATAC2=ATAC_1[11085:(11084+4313),]
p=pheatmap(ATAC1,cluster_col=FALSE)
ATAC1_1=ATAC1[p$tree_row$order,]
RNA1_1=EXP1[p$tree_row$order,]
p=pheatmap(ATAC2,cluster_col=FALSE)
ATAC2_1=ATAC2[p$tree_row$order,]
RNA2_1=EXP2[p$tree_row$order,]

ATAC_f=rbind(ATAC1_1,ATAC2_1)
EXP_f=rbind(RNA1_1,RNA2_1)
color1=colorRampPalette(c("#5ab15f",'#a3d700','white',"#faf513","#f1b319"))(100)

p1=pheatmap(EXP_f,show_rownames=FALSE,cluster_col=FALSE,cluster_row=FALSE,color=color1,treeheight_row=0,treeheight_col=0,gaps_row=c(11084))
pdf('DRUG4.pdf',width=5,height=10)
p1
dev.off()

p=pheatmap(ATAC_f,color=color0,cluster_col=FALSE,cluster_row=FALSE,show_rownames=FALSE,treeheight_row=0,treeheight_col=0,gaps_row=c(11084))
pdf('ATAC4.pdf',width=5,height=10)
p
dev.off()
