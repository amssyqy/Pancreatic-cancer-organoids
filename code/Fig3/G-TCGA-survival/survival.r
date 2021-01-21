library(cgdsr)
# setwd('E:\\mine\\work\\pancreatic\\version3\\code\\Fig3\\G-TCGA-survival')
mycgds = CGDS("https://www.cbioportal.org/")
# Get list of cancer studies at server
# Get available case lists (collection of samples) for a given cancer study
mycancerstudy = getCancerStudies(mycgds)[203,1] #"paad_tcga" 
mycaselist = getCaseLists(mycgds,mycancerstudy)[1,1]
# Get available genetic profiles
# Get data slices for a specified list of genes, genetic profile and case list

# Get clinical data for the case list
myclinicaldata = getClinicalData(mycgds,mycaselist)
## get mutation data  
# Create CGDS object先画图
## get mutation data  paad_tcga_rna_seq_v2_mrna
rankAA = read.table('NMF-result.txt',sep='\t',header=TRUE,row.names=1)

subtype = as.character(1:158)
AA1=t(rankAA)
for (i in 1:length(subtype))
{
	if (i<=17){
	subtype[i] = "Progenitor"
	}else if (i<=55){
	subtype[i] = "Classical"
	}else if (i<=103){
	subtype[i] = "Basal-like"
	}else {
	subtype[i]= "Glycomet"
	}
}

subtype=as.factor(subtype)
index0=c(1:158)
subtype=subtype[index0]
library(survival)
choose_columns=c('OS_STATUS','OS_MONTHS')
choose_clinicaldata=myclinicaldata[,choose_columns]
dat=cbind(choose_clinicaldata[rownames(AA1)[index0],c('OS_STATUS','OS_MONTHS')],subtype)
A=dat$OS_MONTHS > 0
dat=dat[A,]
B=!is.na(dat$OS_STATUS)
dat=dat[B,]
dat$OS_STATUS=as.character(dat$OS_STATUS) 
attach(dat)
library("survminer")
ann_colors = c("#00a087","#4dbbd5","#8491b4","#f39b7f")


my.surv <- Surv(OS_MONTHS,OS_STATUS=="1:DECEASED")
kmfit1 <- survfit(my.surv~subtype,data = dat)
A=ggsurvplot(kmfit1,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE,palette= ann_colors)
pdf("ggsurvplot.pdf",width=6, height=8 ,onefile = FALSE)
A
dev.off()

