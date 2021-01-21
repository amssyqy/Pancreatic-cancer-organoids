col_cut=c(7,18,27)
library(pheatmap)
data0=read.table("finalheatmap.txt",sep='\t',header=TRUE)
anno1=c(rep("Class1",7),rep("Class2",11),rep("Class3",9),rep("Class4",13))
anno1=data.frame(anno1)
rownames(anno1)=as.character(t(colnames(data0))[2:41])
ann_colors = list(
    anno1 = c(Class1 = "#4dbbd5", Class2 = "#00a087",Class3 = "#f39b7f",Class4="#8491b4"), anno=c("Class1"= "#4dbbd5","Class2"= "#00a087","Class3"= "#f39b7f","Class4"="#8491b4"))
H=6
W=5
library(RColorBrewer)
color1=colorRampPalette(c("#3363aa","#1384d5","#23b2ae"))(4)
color10=colorRampPalette(c(color1[1],color1[2]))(35)
color11=colorRampPalette(c(color1[2],color1[4]))(15)
color2=colorRampPalette(c( "#23b2ae","#c5bc5e","#faf513"))(4)
color20=colorRampPalette(c(color2[1],color2[3]))(15)
color21=colorRampPalette(c(color2[3],color2[4]))(35)
mycolor=c(color10,color11,color20,color21)
anno=data.frame(data0$anno)
rownames(data0)=as.character(1:dim(data0)[1])
rownames(anno)=rownames(data0)
colnames(anno)="anno"
pdf(paste0("finalheatmap.pdf"),height=H,width=W)
pheatmap(data0[2:41],scale="row",color = mycolor,fontsize=7,fontsize_row = 70*H/dim(data0)[1], fontsize_col = 7,cluster_cols=FALSE, cluster_rows=FALSE,annotation_col=anno1,gaps_col = col_cut,gaps_row = c(288,755,1099), annotation_row=anno,annotation_colors = ann_colors)
dev.off()
