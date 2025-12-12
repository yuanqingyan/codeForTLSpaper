######################
.libPaths(c("/projects/p32215/tools/forR4.4Seu5.1.0",'/home/yyw9094/R/forR4.4Seu5.1.0',"/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4"))
######################

library("Seurat")#
packageVersion("Seurat") 
library(future)
library(UCell)
library(ggplot2)
library(magrittr) ###

saveFolder<-"/projects/b1146/BharatLab/Yuanqing/Projects/Spatial/IPF"
FigOut<-"/projects/b1042/Yuanqing/Sptial/data/output"
folderPath<-"/projects/b1146/BharatLab/20240530__190352__20240530_X1_Bharat"
(sampleFile<-list.files(folderPath,pattern="__"))
(pdata<-data.frame(sam=sampleFile,disease=c("IPF3R", "IPF5", "IPF4", "Donor2", "IPF3L", "IPF2", "IPF1", "Donor1")))

source("source___studyNiche_allData.R")
##########################################################################################################################################
###################################                                                                 ######################################
###################################                                                                 ######################################
###################################                                                                 ######################################
##########################################################################################################################################

(allF<-list.files(saveFolder,pattern="annotated.RDS"))
IPF0<-readRDS(sprintf("%s/IPF_mergData_annotated.RDS",saveFolder))
IPF<-subset(IPF0, CellType %in% 'unknown', invert=T)



##############check the coordinate of IPF5 ########
head(IPF@images$fov_IPF5$centroids@coords)

set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                            group.by="CellType",
                            assay = "niche6", 
                            cluster.name = "niches6",
                            buffer=500,
                            neighbors.k = 30, 
                            niches.k = 6,
                            Kiter.max=50)

set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                            group.by="CellType",
                            assay = "niche7", 
                            cluster.name = "niches7",
                            buffer=500,
                            neighbors.k = 30, 
                            niches.k = 7,
                            Kiter.max=50)

set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                            group.by="CellType",
                            assay = "niche8", 
                            cluster.name = "niches8",
                            buffer=500,
                            neighbors.k = 30, 
                            niches.k = 8,
                            Kiter.max=50)
set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                            group.by="CellType",
                            assay = "niche9", 
                            cluster.name = "niches9",
                            buffer=500,
                            neighbors.k = 30, 
                            niches.k = 9,
                            Kiter.max=50)

set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                            group.by="CellType",
                            assay = "niche10", 
                            cluster.name = "niches10",
                            buffer=500,
                            neighbors.k = 30, 
                            niches.k = 10,
                            Kiter.max=50)
set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                            group.by="CellType",
                            assay = "niche11", 
                            cluster.name = "niches11",
                            buffer=500,
                            neighbors.k = 30, 
                            niches.k = 11,
                            Kiter.max=50)
set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                            group.by="CellType",
                            assay = "niche12", 
                            cluster.name = "niches12",
                            buffer=500,
                            neighbors.k = 30, 
                            niches.k = 12,
                            Kiter.max=50)
set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                            group.by="CellType",
                            assay = "niche13", 
                            cluster.name = "niches13",
                            buffer=500,
                            neighbors.k = 30, 
                            niches.k = 13,
                            Kiter.max=50)
set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                      group.by="CellType",
                      assay = "niche14", 
                      cluster.name = "niches14",
                      buffer=500,
                      neighbors.k = 30, 
                      niches.k = 14,
                      Kiter.max=50)
set.seed(123)
IPF<-FindNicheAllData(object=IPF, 
                      group.by="CellType",
                      assay = "niche15", 
                      cluster.name = "niches15",
                      buffer=500,
                      neighbors.k = 30, 
                      niches.k = 15,
                      Kiter.max=50)

saveRDS(IPF,file=sprintf("%s/IPF_withNichAllData.RDS",saveFolder))


####################################################################################
##########     outptu niche information       ######################################
####################################################################################
IPF_Nich<-readRDS(file=sprintf("%s/IPF_withNichAllData.RDS",saveFolder))# ###6~15 niches
table(IPF_Nich$niches12,IPF_Nich$orig.ident)
table(IPF_Nich$CellType,IPF_Nich$niches12)

outMeta<-IPF_Nich@meta.data
outMeta$SamCellID<-row.names(outMeta)
outMeta$CellID<-sapply(strsplit(row.names(outMeta),split="\\_"),function(x) x[2]);head(outMeta)
write.csv(outMeta,file=sprintf("%s/XeniumMeta_ILD_B_Paper.csv",saveFolder))

