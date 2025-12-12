######################
.libPaths(c("/projects/p32215/tools/forR4.4Seu5.1.0",'/home/yyw9094/R/forR4.4Seu5.1.0',"/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4"))
######################

library("Seurat")
packageVersion("Seurat") ###this should be 5.1.0
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
source("SpatialSourceCode.R") #
##########################################################################################################################################
###################################                                                                 ######################################
##########################################################################################################################################

IPF<-readRDS(file=sprintf("%s/IPF_mergData_annotated.RDS",saveFolder))
IPF_F<-subset(IPF, CellType %in% 'unknown', invert=T);unique(IPF_F$CellType)
DimPlot(IPF_F,group.by="CellType",label=T,repel=T,raster=T)+NoLegend()
(tabCT_PT<-as.data.frame.matrix(table(IPF_F$CellType, IPF_F$orig.ident)))

##########################################################################################################################################
###################################                                                                 ######################################
##########################################################################################################################################
DefaultAssay(IPF_F)<-'Xenium';Idents(IPF_F)<-'CellType'
CTLevel<-c('AT2','AT1','Ciliated','Club','Basal','KRT17+KRT5-',
           'En',"Capillary.En","Artery.En","Lym.En",
           "Fibro","AlvF","MyoF","SMC","Pericyte","Mesothelial",
           "B","Plasma","CD4T","CD8T","NK",
           "Monocyte","DC","Macro","MoM","IM","Neutrophils","Chemokine.Imm","Prolif.Cells","Mast"
           );length(CTLevel)
IPF_F$CellType<-factor(IPF_F$CellType,levels=rev(CTLevel))
Idents(IPF_F)<-"CellType"

pltMarkerGene<-c("HHIP", "AGER","FOXJ1", "CTSE",'KRT17','KRT5',
                 'VWF','TMEM100','FCN3','GJA5','MMRN1',
                 "PDGFRA","SFRP2",'FGFR4','ACTA2','ASPN','MYH11',"COX4I2",'UPK3B',
                 'CD79A','MS4A1','MZB1',
                 'CD3E','CD40LG',"CD8A","KLRD1",
                 'CD300E','FCN1',
                 "CLEC9A",'S100B',
                 'MARCO','SPP1',"F13A1","FOLR2",
                 "FCGR3B","CSF3R","CCL19","CXCL9",
                 'MKI67','MS4A2');length(pltMarkerGene)

(dotT<-DotPlot(object =IPF_F, 
              features = pltMarkerGene,
              cols = c("lightgrey", "red"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1)))


pdf(sprintf("%s/CellType_Dotplot.pdf",FigOut),width=7.87*1.5,height=7.87, family="ArialMT")
print(dotT)
dev.off()


#########
library(RColorBrewer)
colors_df <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(IPF_F$CellType)))
colors_df<-c("#E41A1C", "#B43547", "#845172", "#546C9D", "#3983AC", "#3F908E", "#FFB315", "#FFD723", "#FFFA31", "#E9D630",
              "#85658D", "#9B4F9D", "#B75D70", "#D46A43", "#F07816", "#FF9007",
              'blue','cyan','green1','green4',####. B-CD8T   "#FFB315", "#FFD723", "#FFFA31", "#E9D630",
             "#D0A72D", "#B8782A", "#AB5832", "#C1645C",
             "#D87085", "#EE7CAF", "#E685B8", "#CC8BAD", "#B292A3", "#999999")

IPF_F$CellType<-factor(IPF_F$CellType,levels=CTLevel)
colSch<-data.frame(coldf=colors_df,ct=CTLevel);head(colSch)
colSch[colSch$ct %in% c("B",'Plasma','CD4T','CD8T'),]

dp<-DimPlot(IPF_F,group.by="CellType",label=F,repel=T,raster=T,cols=colors_df) & NoAxes()
dp1<-DimPlot(IPF_F,group.by="CellType",label=T,repel=T,raster=T,cols=colors_df)

DimPlot(IPF_F, label=T, group.by="CellType", raster=F,
        cells.highlight= list(row.names((subset(IPF_F, CellType == 'Chemokine.Imm')@meta.data))),
        cols.highlight = c("red"), cols= "grey")+NoLegend()

DimPlot(IPF_F, label=T, group.by="CellType", raster=F,
        cells.highlight= list(row.names((subset(IPF_F, CellType == 'Prolif.Cells')@meta.data))),
        cols.highlight = c("green"), cols= "grey")+NoLegend()

pdf(sprintf("%s/CT_Dimplot_labeled.pdf",FigOut),width=7.87*1.5,height=7.87, family="ArialMT")
print(dp1)
dev.off()

tiff(sprintf("%s/CT_Dimplot_labeled.tiff",FigOut),width=7.87*1.5,height=7.87, res=600, units="in",compression='lzw')
print(dp)
dev.off()


#######################
ImageDimPlot(IPF_F, fov='fov_Donor1', group.by = "CellType", size = 0.6, cols=c(rep('grey',0),'red',rep('grey',29)), dark.background = F)+NoLegend() + ggtitle("test")
table(IPF_F$CellType,IPF_F$orig.ident)
ImageDimPlot(IPF_F, fov='fov_Donor2', group.by = "CellType", size = 0.6, cols=c(rep('grey',0),'red',rep('grey',29)), dark.background = F)+NoLegend() + ggtitle("test")
#######################
require(gridExtra)
table(IPF_F$CellType, IPF_F$orig.ident)
IPF_F$CellType<-factor(IPF_F$CellType,levels=CTLevel);unique(IPF_F$CellType)

(uniID<-unique(IPF_F$orig.ident))
lapply(4:length(uniID),function(iD){
  idname<-uniID[iD];print(idname)
  allPlot_id<-lapply(1:length(CTLevel),function(x){
    ct<-CTLevel[x]
    start<-(x-1)
    end<-(30-x)
    dimDonor1<-ImageDimPlot(IPF_F, 
                            fov=sprintf('fov_%s',idname),
                            group.by = "CellType", 
                            size = 0.6,
                            cols=c(rep('grey',start),'red',rep('grey',end)),
                            dark.background = F)+NoLegend()+ggtitle(sprintf("%s",ct))
    return(dimDonor1)
  })
  pdf(sprintf("%s/%s_imgDimplot.pdf",FigOut,idname),width=7.87*3,height=7.87*5, family="ArialMT")
  do.call(grid.arrange, c(allPlot_id, ncol=6,nrow=5))
  dev.off()
})



allPlot<-lapply(1:length(CTLevel),function(x){
  ct<-CTLevel[x]
  start<-(x-1)
  end<-(30-x)
  dimDonor1<-ImageDimPlot(IPF_F, 
                          fov='fov_Donor1',
                          group.by = "CellType", 
                          size = 0.6,
                          cols=c(rep('grey',start),'red',rep('grey',end)),
                          dark.background = F)+NoLegend()#+ggtitle(sprintf("%s",ct))
   tiff(sprintf("%s/CellType/ImgDim_D1_%s.tiff",FigOut, ct),width=7.87/1.3,height=7.87, res=600, units="in",compression='lzw')
   print(dimDonor1);dev.off()
})




tiff(sprintf("%s/Donor1_imgDimplot.tiff",FigOut),width=7.87*3,height=7.87*5, family="ArialMT")
do.call(grid.arrange, c(allPlot, ncol=6,nrow=5))
dev.off()

###donor2
allPlot2<-lapply(1:length(CTLevel),function(x){
  ct<-CTLevel[x]
  start<-(x-1)
  end<-(30-x)
  dimDonor1<-ImageDimPlot(IPF_F, 
                          fov='fov_Donor2',
                          group.by = "CellType", 
                          size = 0.6,
                          cols=c(rep('grey',start),'red',rep('grey',end)),
                          dark.background = F)+NoLegend()#+ggtitle(sprintf("%s",ct))
  tiff(sprintf("%s/CellType/ImgDim_D2_%s.tiff",FigOut, ct),width=7.87/1.3,height=7.87, res=600, units="in",compression='lzw')
  print(dimDonor1);dev.off()
#  return(dimDonor1)
})
# pdf(sprintf("%s/Donor2_imgDimplot.pdf",FigOut),width=7.87*3,height=7.87*4.5, family="ArialMT")
# do.call(grid.arrange, c(allPlot2, ncol=6,nrow=5))
# dev.off()

###############
Dor1BP<-ImageDimPlot(IPF_F, fov='fov_Donor1', group.by = "CellType", size = 0.6, 
                     cols=c(rep('grey',16),colors_df[17:18],rep('grey',12)), dark.background = F)+NoLegend()
tiff(sprintf("%s/CellType/ImgDim_Dor1_BP.tiff",FigOut),width=7.87/1.3,height=7.87, res=600, units="in",compression='lzw')
print(Dor1BP);dev.off()
Dor1T<-ImageDimPlot(IPF_F, fov='fov_Donor1', group.by = "CellType", size = 0.6, 
                     cols=c(rep('grey',18),colors_df[19:20],rep('grey',10)), dark.background = F)+NoLegend()
tiff(sprintf("%s/CellType/ImgDim_Dor1_T.tiff",FigOut),width=7.87/1.3,height=7.87, res=600, units="in",compression='lzw')
print(Dor1T);dev.off()


IPF2BP<-ImageDimPlot(IPF_F, fov='fov_IPF2',
             group.by = "CellType", size = 0.6, cols=c(rep('grey',16),colors_df[17:18],rep('grey',12)), dark.background = F)+NoLegend()
tiff(sprintf("%s/CellType/ImgDim_IPF2_BP.tiff",FigOut),width=7.87/1.3,height=7.87, res=600, units="in",compression='lzw')
print(IPF2BP);dev.off()
IPF2T<-ImageDimPlot(IPF_F, fov='fov_IPF2',
                     group.by = "CellType", size = 0.6, cols=c(rep('grey',18),colors_df[19:20],rep('grey',10)), dark.background = F)+NoLegend()
tiff(sprintf("%s/CellType/ImgDim_IPF2_T.tiff",FigOut),width=7.87/1.3,height=7.87, res=600, units="in",compression='lzw')
print(IPF2T);dev.off()



IPF4BP<-ImageDimPlot(IPF_F, fov='fov_IPF4',
                    group.by = "CellType", size = 0.6, cols=c(rep('grey',16),colors_df[17:18],rep('grey',12)), dark.background = F)+NoLegend()
tiff(sprintf("%s/CellType/ImgDim_IPF4_BP.tiff",FigOut),width=7.87/1.3,height=7.87, res=600, units="in",compression='lzw')
print(IPF4BP);dev.off()

IPF4T<-ImageDimPlot(IPF_F, fov='fov_IPF4',
                     group.by = "CellType", size = 0.6, cols=c(rep('grey',18),colors_df[19:20],rep('grey',10)), dark.background = F)+NoLegend()
tiff(sprintf("%s/CellType/ImgDim_IPF4_T.tiff",FigOut),width=7.87/1.3,height=7.87, res=600, units="in",compression='lzw')
print(IPF4T);dev.off()
