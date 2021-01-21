##figureS2A draw with Integrative Genomics Viewer(IGV)
##figureS2B draw with GraphPad Prism
##figureS2C
library(maftools)##2.2.10
library(RColorBrewer)##1.1-2
load("figureS2C_data.RData")
sample_order<-c("CAS-DAC-16-O","CAS-DAC-16-T","CAS-DAC-20-O","CAS-DAC-20-T","CAS-DAC-26-O","CAS-DAC-26-T","CAS-DAC-28-O","CAS-DAC-28-T","CAS-DAC-29-O","CAS-DAC-29-T","CAS-DAC-33-O","CAS-DAC-33-T","CAS-DAC-34-O","CAS-DAC-34-T","CAS-DAC-42-O","CAS-DAC-42-T","CAS-NEN-1-O","CAS-NEN-1-T","CAS-ACC-1-O","CAS-ACC-1-T")
col = brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins','In_Frame_Ins', 'In_Frame_Del',"Splice_Site")
##draw oncoplot
oncoplot(maf =figureS2C_data, colors = col, drawRowBar = F,drawColBar = F,
         genes = c("KRAS","TP53","SMAD4","CDKN2A","RNF43","ARID1A"),
         sortByMutation   = TRUE, annotationFontSize = 1,
         legendFontSize = 1,removeNonMutated = F,showTumorSampleBarcodes =T 
         ,sampleOrder = sample_order,barcode_mar = 7,showTitle = F)
##figureS2D
library(RCircos)##1.2.1
library(Cairo)##1.5-11
data(UCSC.HG19.Human.CytoBandIdeogram)
load("CNV_tissue_figureS2D.RData")
load("CNV_organoid_figureS2D.RData")
for (sample in colnames(CNV_tissue_figureS2D)[1:5])
{ 
  CNA_draw_organoid<-CNV_organoid_figureS2D[,c("Chromosome","chromStart","chromEnd","GeneName",sample)]
  CNA_draw_organoid[,5]<-as.numeric(CNA_draw_organoid[,5])
  CNA_draw_organoid_amp<-CNA_draw_organoid[as.numeric(CNA_draw_organoid[,5])>=0,]
  CNA_draw_organoid_del<-CNA_draw_organoid[as.numeric(CNA_draw_organoid[,5])<=0,]
  CNA_draw_tissue<-CNV_tissue_figureS2D[,c("Chromosome","chromStart","chromEnd","GeneName",sample)]
  CNA_draw_tissue[,5]<-as.numeric(CNA_draw_tissue[,5])
  CNA_draw_tissue_amp<-CNA_draw_tissue[as.numeric(CNA_draw_tissue[,5])>=0,]
  CNA_draw_tissue_del<-CNA_draw_tissue[as.numeric(CNA_draw_tissue[,5])<=0,]
  tiff(file=paste0(sample,"_circos.tiff"),compression="lzw",units="in",res=1000,pointsize=8,height=4.5,width=4.5)
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
  RCircos.Set.Core.Components(cyto.info, chr.exclude=c("chrX","chrY"),tracks.inside=4, tracks.outside=0 )
  RCircos.Set.Plot.Area()     
  RCircos.Chromosome.Ideogram.Plot() 
  param<-RCircos.Get.Plot.Parameters()
  param$track.background<-"ghostwhite"
  param$point.size=1
  param$grid.line.color=NULL
  param$track.height=0.16
  param$track.padding=0.04
  param$heatmap.width=10
  RCircos.Reset.Plot.Parameters(param)
  ##CNA heatmap
  RCircos.Heatmap.Plot(CNA_draw_organoid_amp, 5, 1, "in",min.value = -2,max.value = 2)
  RCircos.Heatmap.Plot(CNA_draw_organoid_del, 5, 3, "in",min.value = -2,max.value = 2)
  RCircos.Heatmap.Plot(CNA_draw_tissue_amp, 5, 2, "in",min.value = -2,max.value = 2)
  RCircos.Heatmap.Plot(CNA_draw_tissue_del, 5, 4, "in",min.value = -2,max.value = 2)
  dev.off()
}
