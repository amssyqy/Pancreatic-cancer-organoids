col_cut=c(7,18,27)
library(pheatmap)
data0=read.table("gly.txt",sep='\t',header=FALSE,row.names=1)
subtype=c(rep("Classical-like",7),rep("basal-like",11),rep("Classical-Progenitor",9),rep("Class4",13))
subtype=data.frame(subtype)
rownames(subtype)=as.character(t(colnames(data0))[1:40])
ann_colors = list(subtype = c("Classical-like" = "#4dbbd5", "basal-like" = "#00a087","Classical-Progenitor" = "#f39b7f",Class4="#8491b4"))
H=2
W=6
library(RColorBrewer)
ATAC_zscore=as.matrix(data0)
data0=as.matrix(data0)
for (i in 1:dim(data0)[1])
{
  mean0=mean(data0[i,])
  sigma=sd(data0[i,])
  ATAC_zscore[i,]=(data0[i,]-mean0)/sigma
}
ATAC_zscore[ATAC_zscore>2]=2
color1=colorRampPalette(c("#3363aa","#1384d5","#23b2ae","#c5bc5e","#faf513"))(100)

pdf(paste0("gly1.pdf"),height=H,width=W)
pheatmap(ATAC_zscore,color=color1,fontsize=7,treeheight_row=10, border_color =NA,fontsize_row = 70*H/dim(data0)[1], fontsize_col = 0.1,cluster_cols=FALSE,annotation_col=subtype,gaps_col = col_cut, annotation_colors = ann_colors)
dev.off()