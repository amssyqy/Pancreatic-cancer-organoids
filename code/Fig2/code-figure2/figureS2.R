##figureS2a draw with Integrative Genomics Viewer(IGV)
##figureS2b/g draw with GraphPad Prism
##figureS2c/d
library(MutationalPatterns)
library(BSgenome)
library(CancerSubtypes)
library(ComplexHeatmap)
library(RColorBrewer)
library(gridExtra)
library(NMF)
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
library(ref_genome, character.only = TRUE)
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome, type = "all")
snv <- get_mut_type(vcfs, type = "snv")
mut_mat <- mut_matrix(vcf_list = snv, ref_genome = ref_genome)
mut_mat <- mut_mat + 0.0001
estimate <- nmf(mut_mat, rank = 2:8, method = "brunet", 
                nrun = 1000, seed = 123456, .opt = "v-p")
nmf_res <- extract_signatures(mut_mat, rank = 3, nrun = 1000, single_core = TRUE)
colnames(nmf_res$signatures) <- c("Signature A", "Signature B","Signature C")
rownames(nmf_res$contribution) <- c("Signature A", "Signature B","Signature C")
plot_96_profile(nmf_res$signatures,condensed = T,colors = RColorBrewer::brewer.pal(6,"Set2"))
########cancer_signatures,cosmic v2
sort(apply(cancer_signatures,2,function(x) cos_sim(nmf_res$signatures[,1],x)),decreasing = T)
sort(apply(cancer_signatures,2,function(x) cos_sim(nmf_res$signatures[,2],x)),decreasing = T)
sort(apply(cancer_signatures,2,function(x) cos_sim(nmf_res$signatures[,3],x)),decreasing = T)
########distribution
plot_contribution(nmf_res$contribution, nmf_res$signature, mode = "relative",palette = RColorBrewer::brewer.pal(8,"Pastel1")[2:4])
#######signature_clustering
nmf_matrix<-nmf_res$contribution
nmf_matrix<-t(nmf_matrix)
result=ExecuteCC(clusterNum=3,d=as.matrix(t(nmf_matrix)),maxK=6,clusterAlg="hc",distance="pearson",title="hcpcc-3",innerLinkage="ward.D2",reps=5000)
#load("data-figure2/figureS2d_data.RData")#
sample_anno<-data.frame(colnames(t(nmf_matrix)),group=result$group,row.names = 1)
annotation_color<-RColorBrewer::brewer.pal(8,"Pastel1")[2:4]
names(annotation_color)=unique(sample_anno$group)
ComplexHeatmap::pheatmap(t(nmf_matrix),
                            name="Weight",scale = "column",
                            col = colorRampPalette(c("royalblue","white","red"))(30000),
                            cluster_cols = T, cellwidth = 7,cluster_rows = F,
                            cellheight = 7.5,border_color = "white",
                            clustering_distance_rows = "pearson",
                            show_colnames = T,annotation_col = sample_anno,
                            show_rownames =  T,annotation_colors  = list(group=annotation_color),
                            column_split = factor(sample_anno$group,levels = c(1:4)),
                            show_row_dend = F,
                            fontsize = 8)

##figureS2e/f
library(VennDiagram)
venn.diagram(x=list(Tissue, Organoid)
             ,filename ="sample_venn.png"
             ,height=2400
             ,width= 2400
             ,resolution =600
             ,imagetype="png"
             ,col="black"
             ,fill=c("cornflowerblue", "darkorange1")
             ,alpha=c(0.6, 0.6)
             ,lwd=c(0.2, 0.2)
             ,cex=0.12
             ,print.mode = "raw"
             ,main = sampleID
             ,main.cex = 0.13)

##figureS2h
library(maftools)##2.2.10
library(RColorBrewer)##1.1-2
load("data-figure2/figureS2h_data.RData")
sample_order<-c("CAS-DAC-16-O","CAS-DAC-16-T","CAS-DAC-20-O","CAS-DAC-20-T","CAS-DAC-26-O","CAS-DAC-26-T","CAS-DAC-28-O","CAS-DAC-28-T","CAS-DAC-29-O","CAS-DAC-29-T","CAS-DAC-33-O","CAS-DAC-33-T","CAS-DAC-34-O","CAS-DAC-34-T","CAS-DAC-42-O","CAS-DAC-42-T","CAS-NEN-1-O","CAS-NEN-1-T","CAS-ACC-1-O","CAS-ACC-1-T")
col = brewer.pal(n = 8, name = 'Paired')
names(col) = c('Frame_Shift_Del','Missense_Mutation', 'Nonsense_Mutation', 'Multi_Hit', 'Frame_Shift_Ins','In_Frame_Ins', 'In_Frame_Del',"Splice_Site")
##draw oncoplot
oncoplot(maf =figureS2h_data, colors = col, drawRowBar = F,drawColBar = F,
         genes = c("KRAS","TP53","SMAD4","CDKN2A","RNF43","ARID1A"),
         sortByMutation   = TRUE, annotationFontSize = 1,
         legendFontSize = 1,removeNonMutated = F,showTumorSampleBarcodes =T 
         ,sampleOrder = sample_order,barcode_mar = 7,showTitle = F)

##figureS2i
library(RCircos)##1.2.1
library(Cairo)##1.5-11
data(UCSC.HG19.Human.CytoBandIdeogram)
load("data-figure2/figureS2i_data.RData")
for (sample in colnames(CNV_tissue_figureS2i)[1:5])
{ 
  CNA_draw_organoid<-CNV_organoid_figureS2i[,c("Chromosome","chromStart","chromEnd","GeneName",sample)]
  CNA_draw_organoid[,5]<-as.numeric(CNA_draw_organoid[,5])
  CNA_draw_organoid_amp<-CNA_draw_organoid[as.numeric(CNA_draw_organoid[,5])>=0,]
  CNA_draw_organoid_del<-CNA_draw_organoid[as.numeric(CNA_draw_organoid[,5])<=0,]
  CNA_draw_tissue<-CNV_tissue_figureS2i[,c("Chromosome","chromStart","chromEnd","GeneName",sample)]
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
