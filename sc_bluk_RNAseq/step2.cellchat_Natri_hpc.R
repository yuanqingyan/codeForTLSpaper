######################
.libPaths(c("/projects/p32215/tools/forSlurmR440","/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4","/projects/p32215/tools/forR4.4Seu5.1.0"))
######################
library(Seurat)
library(harmony)
library(cowplot)
library(dplyr)
library(ggplot2)

FigOut<-"/projects/b1042/Yuanqing/SC/ILD/output"


list.files("./PubD/GSE227136")
sc<-readRDS("./PubD/GSE227136/GSE227136_ILD_all_celltypes_Seurat.rds");dim(sc)
DimPlot(sc,group.by='manual_annotation_1',label=T,repel=T,raster=T)+NoLegend()
bsub<-readRDS("./PubD/GSE227136/RDS/Natri.B_Subtype.RDS")
head(bsub)

sc$CellType<-sc$manual_annotation_1
test<-data.frame(test1=sapply(strsplit(row.names(sc@meta.data),split="\\_"),function(x) x[3]),
                 test2=sc$lineage);table(test)
sc$CellID<-gsub("_4","",row.names(sc@meta.data))

metaF<-bsub@meta.data
uniqFCT<-unique(metaF$FineCT);uniqFCT
for(i in 1:length(uniqFCT)){
  tempCT<-uniqFCT[i]
  temp<-metaF[metaF$FineCT %in% tempCT,]
  sc$CellType<-ifelse(sc$CellID %in% row.names(temp),tempCT,sc$CellType)
}
saveRDS(sc,file=sprintf("%s/sc_____clean.RDS",FigOut))

sc<-readRDS(file=sprintf("%s/sc_____clean.RDS",FigOut))
DimPlot(sc,group.by='CellType',label=T,repel=T,raster=T)+NoLegend()
DimPlot(subset(sc,manual_annotation_1 %in% c("B cells",'Plasma')),
        group.by='CellType',
        label=T,
        repel=T,
        raster=F)+NoLegend()
table(sc$CellType)
table(sc$manual_annotation_1)


# #########################################################
library(dplyr)
library(Seurat)
library(CellChat)
library(patchwork)
library(singleGEO)
library(reshape2)
library(BioCircos)
library(RColorBrewer)
require("gplots")
source("./sourceCellChat.R")


unique(sc$CellType)
################################
seu_Dis<-subset(sc,Status %in% "Control",invert=TRUE);unique(seu_Dis$Status)
ctTable<-table(seu_Dis$CellType)
selCT<-ctTable[ctTable>100];length(selCT)

seuObjAll<-subset(seu_Dis,CellType %in% names(selCT))
saveRDS(seuObjAll,file=sprintf("%s/seuObjAll.RDS",FigOut))

CC_Natri_All<-CellChat_pipeline(sub_CC=seuObjAll,type="truncatedMean",trim=0.25)
saveRDS(CC_Natri_All,file="./PubD/GSE227136/RDS/CellChat_Natri_AllCellType.RDS")

