col_cut=c(38,100,123)
library(pheatmap)
# setwd("E:\\mine\\work\\pancreatic\\version3\\code\\Fig3\\F-heatmap-TCGA")
data0=read.table("TCGA-overlap-heatmap.txt",sep='\t',header=FALSE)
anno1=c(rep("class1",38),rep("class2",62),rep("class3",23),rep("class4",56))
anno1=data.frame(anno1)
rownames(anno1)=as.character(t(colnames(data0))[1:179])
anno=c(rep('class1',106),rep('class2',81),rep('class3',94),rep('class4',78))
anno=data.frame(anno)
rownames(anno)=1:359
ann_colors = list(anno1 = c(class1 = "#4dbbd5", class2 = "#00a087",class3 = "#f39b7f",class4="#8491b4"), anno=c(class1 = "#4dbbd5", class2 = "#00a087",class3 = "#f39b7f",class4="#8491b4"))
H=3
W=8
library(RColorBrewer)
color1=colorRampPalette(c("navy", "white"))(6)
color10=colorRampPalette(c(color1[1],color1[2]))(10)
color11=colorRampPalette(c(color1[2],color1[4]))(20)
color12=colorRampPalette(c(color1[4],color1[6]))(20)
color2=colorRampPalette(c( "white", "firebrick3"))(4)
color20=colorRampPalette(color2[1:3])(20)
color21=colorRampPalette(color2[3:4])(30)
mycolor=c(color10,color11,color12,color20,color21)

library(RColorBrewer)
# color1=colorRampPalette(c("navy", "white"))(4)
# color10=colorRampPalette(c(color1[1],color1[2]))(35)
# color11=colorRampPalette(c(color1[2],color1[4]))(15)
# color2=colorRampPalette(c( "white", "firebrick3"))(4)
# color20=colorRampPalette(c(color2[1],color2[3]))(15)
# color21=colorRampPalette(c(color2[3],color2[4]))(35)
# mycolor=c(color10,color11,color20,color21)
# rownames(data0)=as.character(1:dim(data0)[1])
ATAC_zscore=as.matrix(data0)
data0=as.matrix(data0)
for (i in 1:dim(data0)[1])
{
  mean0=mean(data0[i,])
  sigma=sd(data0[i,])
  ATAC_zscore[i,]=(data0[i,]-mean0)/sigma
}
ATAC_zscore[ATAC_zscore>0]=sqrt(ATAC_zscore[ATAC_zscore>0])
rownames(ATAC_zscore)=1:359

pdf(paste0("TCGA-overlap-heatmap.pdf"),height=H,width=W)
pheatmap(ATAC_zscore,fontsize=7,color = mycolor, fontsize_row = 0.1, fontsize_col = 0.1,cluster_cols=FALSE,cluster_rows=FALSE,annotation_col=anno1,annotation_row=anno, annotation_colors = ann_colors,gaps_col = col_cut,gaps_row = c(106,187,281))
dev.off()


