######################
.libPaths(c("/projects/p32215/tools/forSlurmR440","/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4","/projects/p32215/tools/forR4.4Seu5.1.0"))
######################
library(Seurat)
library(harmony)
library(cowplot)
library(dplyr)
library(ggplot2)
library(CellChat)
library(patchwork)
library(singleGEO)
library(reshape2)
library(BioCircos)
library(RColorBrewer)
require("gplots")
source("./sourceCellChat.R")

FigOut<-"/projects/b1042/Yuanqing/SC/ILD/output"
###########################
CC_Natri_All<-readRDS(file="./PubD/GSE227136/RDS/CellChat_Natri_AllCellType.RDS")
head(CC_Natri_All$df.net)
unique(CC_Natri_All$df.net$target)
unique(CC_Natri_All$df.net[CC_Natri_All$df.net$source %in% c("Naive.B","Memory.Act.B") &
                          CC_Natri_All$df.net$target %in% c("Adventitial FB",'Alveolar FB',
                                                         'MyoFB','MyoFB - Activated','Pericyte'),
                        c('source','target','interaction_name_2')])

unique(CC_Natri_All$df.net[CC_Natri_All$df.net$source %in% c("Naive.B","Memory.Act.B") &
                             CC_Natri_All$df.net$ligand %in% c("SELL"),
                           c('source','target','interaction_name_2')])

sc<-readRDS(file=sprintf("%s/sc_____clean.RDS",FigOut))
seu_Dis<-subset(sc,Status %in% "Control",invert=TRUE);unique(seu_Dis$Status)
ctTable<-table(seu_Dis$CellType)
selCT<-ctTable[ctTable>100];length(selCT)
seuObjAll<-subset(seu_Dis,CellType %in% names(selCT))

DefaultAssay(seuObjAll)<-"RNA";Idents(seuObjAll)<-"CellType"
VlnPlot(seuObjAll,features=c("BTLA",'GRN','LTA','TNF','MIF','TGFB1'),pt.size=0)##
VlnPlot(seuObjAll,features=c("SELL"),pt.size=0)+NoLegend()

DefaultAssay(sc)<-"RNA";Idents(sc)<-"CellType"
VlnPlot(sc,features=c("SELL"),pt.size=0,split.by='Status')+NoLegend()


##############################################################################################################
#############                 generate the plots  ---panel                                          ##########
##############################################################################################################
cellchat_plot(obj.cellchat=CC_Natri_All$cellchat,
              outFolder=sprintf("%s/NatriCC_all",FigOut),
              objName='NatriCC_all')

(groupSize=as.numeric(table(CC_Natri_All$cellchat@idents)))
netVisual_circle(CC_Natri_All$cellchat@net$count, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(CC_Natri_All$cellchat@net$weight, 
                 vertex.weight = groupSize, 
                 weight.scale = T, 
                 label.edge= F, 
                 title.name = "Interaction weights/strength")

##############################################################################################################
(uniquePathway<-unique(CC_Natri_All$cellchat@netP$pathways));length(uniquePathway) 

head(CC_Natri_All$df.net)
(pltPaW<-unique(CC_Natri_All$df.net[CC_Natri_All$df.net$source %in% c("Naive.B","Memory.Act.B") |
                             CC_Natri_All$df.net$target %in% c("Naive.B","Memory.Act.B"),c('pathway_name')]))


obj.cellchat<-CC_Natri_All$cellchat



outFolder<-FigOut
objName<-"Natri"
 

########################################################################################################
##########################           plot all data for reference             ###########################
########################################################################################################
for(ipath in 1:length(pltPaW)){
  pathways.show <-pltPaW[ipath]
  CC_Natri_All$df.net[CC_Natri_All$df.net$pathway_name==pathways.show,c('source','target','interaction_name','pathway_name')]
  
  for(iloop in 1:2){
    if(iloop==1){tiff(sprintf("%s/%s_%s.VlnExp.tiff",outFolder,objName,pathways.show),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
    if(iloop==2){pdf(sprintf("%s/%s_%s.VlnExp.pdf",outFolder,objName,pathways.show),width=7.1,height=7.1)}
    print(plotGeneExpression(obj.cellchat, signaling = pathways.show, enriched.only = TRUE, type = "violin"))
    dev.off()}
  
  for(iloop in 1:2){
    ##circle plot
    if(iloop==1){tiff(sprintf("%s/%s_%s.circle.tiff",outFolder,objName,pathways.show),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
    if(iloop==2){pdf(sprintf("%s/%s_%s.circle.pdf",outFolder,objName,pathways.show),width=7.1,height=7.1)}
    netVisual_aggregate(obj.cellchat, signaling = pathways.show, layout = "circle")
    dev.off()}


  for(iloop in 1:2){
    ##contribution plot
    if(iloop==1){tiff(sprintf("%s/%s_%s.contribution.tiff",outFolder,objName,pathways.show),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
    if(iloop==2){pdf(sprintf("%s/%s_%s.contribution.pdf",outFolder,objName,pathways.show),width=7.1,height=7.1)}
    conp<-netAnalysis_contribution(obj.cellchat, signaling = pathways.show)
    print(conp)
    dev.off()}
}
  
  
########################################################################################################
##########################           plot selected data for final figure     ###########################
########################################################################################################

targetPathway<-c("CD48",'CLEC','IL16','SELL','BTLA','BAG','BAFF','APRIL','CD23','CD40','CD45','CD46','CXCL','MHC-I','MHC-II','SELPLG')

names(obj.cellchat@net)
unique(obj.cellchat@meta$group)
(unique(obj.cellchat@meta$group))
cellGrp<-c("AT2", "Transitional AT2","AT1", "KRT5-/KRT17+","Secretory - SCGB3A2+","Secretory - SCGB1A1+/SCGB3A2+",
           "Secretory - SCGB1A1+/MUC5B+","Ciliated","Differentiating ciliated","Basal","PNEC","Proliferating - Epi",
           "Arteriole","Systemic venous","Venule","aCap","gCap","Lymphatic",
           "Adventitial FB","Alveolar FB","PLIN2+ FB","Pericyte","MyoFB","MyoFB - Activated","SMC","Mesothelial",
           "Inflamed",
           "Monocyte","Inflammatory monocyte","cDC1","cDC2","moDC","pDC",
           "Monocyte-derived macrophage","Interstitial macrophage","Alveolar macrophage",
           "Proliferating - Imm","Mast",
           "NK","CD8/NKT","CD4","Plasma","B cells","Memory.Act.B","Naive.B")

obj.cellchat2<-obj.cellchat
obj.cellchat2@net$prob <- obj.cellchat2@net$prob[cellGrp, cellGrp, ]
obj.cellchat2@net$pval <- obj.cellchat2@net$pval[cellGrp, cellGrp, ]
dimnames(obj.cellchat2@net$prob)[[1]]

lapply(1:length(targetPathway),function(ix){
  temp<-targetPathway[ix]
  templot<-netVisual_aggregate(obj.cellchat2, 
                               signaling = temp, 
                               layout = "circle")
  
  pdf(sprintf("%s/SelectPath_%s.circle.pdf",outFolder,temp),width=7.1,height=7.1)
  print(templot)
  dev.off()
})





##################################################################################################
################################# stacked vlnplot of genes #######################################
##################################################################################################


############
library(Seurat)
library(patchwork)
library(ggplot2)

## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj, 
                          feature, 
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  + 
    xlab("") + ylab(feature) + ggtitle("") + 
    theme(legend.position = "none", 
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.title.y = element_text(size = rel(1), angle = 0), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = plot.margin ) 
  return(p)
}

## extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}


## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # <- rotate x-axis labels
          axis.ticks.x = element_line())
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


##################################################################################################
################################# stacked vlnplot of genes #######################################
##################################################################################################

feature1<- c("SELL", "CD34", "PODXL")
feature2<- c("IL16", "CD4")
feature3<- c("TNFSF13","TNFRSF17",'TNFRSF13B')

seuObjAll$CellType<-factor(seuObjAll$CellType,levels=cellGrp)
DefaultAssay(seuObjAll)<-"RNA";Idents(seuObjAll)<-"CellType"
sp1<-StackedVlnPlot(obj = seuObjAll, features = feature1)

pdf(sprintf("%s/SelectPath_SELL.StackVln1.pdf",outFolder),width=7.1,height=7.1*1.5)
print(sp1)
dev.off()


pdf(sprintf("%s/SelectPath_SELL.StackVln2.pdf",outFolder),width=7.1,height=7.1*1.1)
print(StackedVlnPlot(obj = seuObjAll, features = feature2))
dev.off()


pdf(sprintf("%s/SelectPath_SELL.StackVln3.pdf",outFolder),width=7.1,height=7.1*1.5)
print(StackedVlnPlot(obj = seuObjAll, features = feature3))
dev.off()

