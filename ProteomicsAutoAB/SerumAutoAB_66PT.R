
######################
.libPaths(c("/projects/p32215/tools/forSlurmR440","/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4","/projects/p32215/tools/forR4.4Seu5.1.0"))
######################

FigOut<-"/projects/b1042/Yuanqing/manu/temp"
autoAB<-readRDS(file="./autoAB_PF.RDS")

############
datFinal_PF<-autoAB$datFinal_PF
pPlot<-autoAB$pPlot

library(ComplexHeatmap)
myCol<-c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494")
ann_colors = list(
  Disease=c("IPF" = myCol[1],
            "ILD_Not_IPF" =myCol[2]))

top_annotation <- HeatmapAnnotation(
  Disease = factor(pPlot$Disease_two),
  col = ann_colors
)

pltData<-t(datFinal_PF)
colnames(pltData)<-sapply(strsplit(row.names(datFinal_PF),split="\\."),function(x) x[2])

gene_order <- rownames(pltData)
gene_sums <- colSums(datFinal_PF)[gene_order]
right_annotation <- rowAnnotation(
  Count = anno_barplot(gene_sums, border = FALSE, gp = gpar(fill = "blue"))
)

(heatFina<-ComplexHeatmap::Heatmap(
  pltData,
  col = c('grey','deeppink'),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  top_annotation= top_annotation,  # Add the custom top annotation here
  #right_annotation = right_annotation,  # Add row annotations here
  show_column_names = TRUE,
  row_order = NULL,
  column_order = NULL,
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.rect(x, y, width, height, gp = gpar(col = "white", lwd = 0.5, fill = NA))
  }
))

plotName="Heatmap_Autoanibody"
pdf(sprintf("%s/%s.pdf",FigOut, plotName),width=3.5*3.5,height=3.5/2.2, family="ArialMT")
print(heatFina)
dev.off()



