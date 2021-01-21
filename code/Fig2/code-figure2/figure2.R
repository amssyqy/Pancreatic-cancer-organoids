##figure 2 left
##preparing data and annoatation
load("mut_table_left.RData") ##load mutation data
load("figure2_annotation.RData") ##load clinical data
load("anno_color.RData")##load annotation colour
col= c("HOMDEL" = "royalblue", "AMP" = "red", "MUT" = "darkgreen","Noncoding"="goldenrod1")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(0, "pt"), 
              gp = gpar(fill = "gray94", col = "white"))},
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(0, "pt"), 
              gp = gpar(fill = col["HOMDEL"], col = "white"))},
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0, "pt"), h-unit(0, "pt"), 
              gp = gpar(fill = col["AMP"], col ="white"))},
  MUT = function(x, y, w, h) {
    grid.points(x, y,pch=16,size=unit(0.6,"char"),
                gp = gpar(fill = col["MUT"], col = "darkgreen"))},
  Noncoding = function(x, y, w, h) {
    grid.rect(x, y,w-unit(0, "pt"), h-unit(0, "pt"),
              gp = gpar(fill = col["Noncoding"], col = "white"))})
heatmap_legend_param = list(title = "Alternations", at = c("HOMDEL","AMP","MUT","Noncoding"), 
                            labels = c("Deep del","Amplification","Protein-coding","Noncoding-site"))
###print
library(ComplexHeatmap)
oncoPrint(mut_table_left,
          show_column_names = T,
          remove_empty_rows = F,
          row_order = rownames(mut_table_left),
          column_order = colnames(mut_table_left),
          alter_fun = alter_fun, 
          col = col, 
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10), 
          pct_gp = gpar(fontsize = 10),
          heatmap_legend_param = heatmap_legend_param,
          row_split = factor(c(rep("exonic", 13),rep("UTR", 2),rep("up/downstream", 14), rep("ncRNA", 2),rep("intronic",7),rep("intergenic",5),rep("CNV",9)),levels = c("exonic","UTR","up/downstream","ncRNA","intronic","intergenic","CNV")),
          right_annotation = NULL,
          top_annotation = HeatmapAnnotation(
              Histologic.Subtype = figure2_annotation[colnames(mut_table_left),]$Histologic.Subtype,
              Gender = figure2_annotation[colnames(mut_table_left),]$Gender,
              Tumor.Stage = figure2_annotation[colnames(mut_table_left),]$Tumor.Stage,
              Differentiation.Grade = figure2_annotation[colnames(mut_table_left),]$Differentiation.Grade,
              col=anno_color),
          alter_fun_is_vectorized =FALSE)

##figure 2 right
load("mut_table_right.RData") ##load mutation data

###print
library(ComplexHeatmap)
oncoPrint(mut_table_right,
          show_column_names = T,
          remove_empty_rows = F,
          row_order = rownames(mut_table_right),
          column_order = colnames(mut_table_right),
          alter_fun = alter_fun, 
          col = col, 
          row_names_gp = gpar(fontsize = 10),
          column_names_gp = gpar(fontsize = 10), 
          pct_gp = gpar(fontsize = 10),
          heatmap_legend_param = heatmap_legend_param,
          row_split = factor(c(rep("exonic", 13),rep("UTR", 1),rep("intronic",3),rep("intergenic",8),rep("CNV",5)),levels = c("exonic","UTR","intronic","intergenic","CNV")),
          right_annotation = NULL,
          top_annotation = HeatmapAnnotation(
            Histologic.Subtype = figure2_annotation[colnames(mut_table_right),]$Histologic.Subtype,
            Gender = figure2_annotation[colnames(mut_table_right),]$Gender,
            Tumor.Stage = figure2_annotation[colnames(mut_table_right),]$Tumor.Stage,
            WHO.Grade.for.NEN = figure2_annotation[colnames(mut_table_right),]$WHO.Grade.for.NEN,
            col=anno_color),
          alter_fun_is_vectorized =FALSE)
