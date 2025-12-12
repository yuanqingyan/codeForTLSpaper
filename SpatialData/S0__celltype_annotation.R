
library("SeuratObject", lib.loc = "/projects/p32215/tools/AddiRLib")
library("Seurat", lib.loc =  "/projects/p32215/tools/AddiRLib")
packageVersion("Seurat") ###version 5.0.1
library(future)
library(UCell)
library(ggplot2)
library(magrittr) ##
library(RColorBrewer)
library(plotly)
source("SpatialSourceCode.R")


saveFolder<-"/projects/b1146/BharatLab/Yuanqing/Projects/Spatial/IPF"
FigOut<-"/projects/b1042/Yuanqing/Sptial/data/output"
folderPath<-"/projects/b1146/BharatLab/20240530__190352__20240530_X1_Bharat"
(sampleFile<-list.files(folderPath,pattern="__"))
(pdata<-data.frame(sam=sampleFile,disease=c("IPF3R", "IPF5", "IPF4", "Donor2", "IPF3L", "IPF2", "IPF1", "Donor1")))


lungGenePanel<-read.csv("Xenium_hLung_v1_metadata.csv")
lungGenePanel$CellName<-gsub(" ","_",lungGenePanel$Annotation)
final100<-read.csv("additional100Genes.csv",row.names=1)
addGene<-as.data.frame(matrix(1,nrow=100,ncol=ncol(lungGenePanel)))
colnames(addGene)<-colnames(lungGenePanel)
addGene$Gene<-final100$x
finalGene<-rbind(lungGenePanel,addGene);dim(finalGene)

indData<-lapply(1:nrow(pdata),function(isam){
  fileName<-pdata[isam,"sam"];plotN<-pdata[isam,"disease"]
  dataPath<-sprintf("%s/%s",folderPath,fileName)
  xenium.obj <- LoadXenium(dataPath, fov = sprintf("fov_%s",plotN))
  xenium.obj$orig.ident<-plotN
  xenium.obj <- subset(xenium.obj, subset = nCount_Xenium > 2)
  Idents(xenium.obj)<-xenium.obj$orig.ident
  saveRDS(xenium.obj,file=sprintf("%s/tempData/%s_xenium.obj_readData2.RDS",FigOut,plotN))
  return(xenium.obj)
});names(indData)<-pdata$disease
saveRDS(indData,file=sprintf("%s/tempData/IPF_indData2.RDS",FigOut))

mergData<-merge(x = indData[[1]], 
                y = indData[2:length(indData)],
                add.cell.ids = names(indData))
mergData <- SCTransform(mergData, assay = "Xenium")
DefaultAssay(mergData)

saveRDS(mergData,file=sprintf("%s/tempData/IPF_mergData2.RDS",FigOut))
#mergData<-readRDS(file=sprintf("%s/tempData/IPF_mergData2.RDS",FigOut))

mergData <- RunPCA(mergData, npcs = 50, features = rownames(mergData))
mergData <- FindNeighbors(mergData, dims = 1:50)
mergData <- FindClusters(mergData, verbose = FALSE, resolution = 1)
mergData <- RunUMAP(mergData, dims = 1:50)
DefaultAssay(mergData) <- "SCT"
mergData <- FindClusters(mergData, verbose = FALSE, resolution = 0.5)
mergData[['Xenium']]$counts<-mergData@assays$SCT@counts
mergData[['Xenium']]$data<-mergData@assays$SCT@data
mergData$Disease<-ifelse(mergData$orig.ident %in% c('Donor1','Donor2'),"Control","PF")

saveRDS(mergData,file=sprintf("%s/tempData/IPF_mergData_cluster2.RDS",FigOut))
#mergData<-readRDS(file=sprintf("%s/tempData/IPF_mergData_cluster2.RDS",FigOut))
DefaultAssay(mergData)<-'Xenium'


########################################################################################################################
#######################################      for cell type  annotation  ################################################
########################################################################################################################
cluster.marker<-FindAllMarkers(mergData, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25,
                               test.use="wilcox")
saveRDS(cluster.marker,file=sprintf("%s/mergData_IPF_cluster.marker2.RDS",FigOut))
#cluster.marker<-readRDS(file=sprintf("%s/mergData_IPF_cluster.marker2.RDS",FigOut))
colnames(finalGene)[1]<-'gene'

library(dplyr)
clusterCell<-dplyr::left_join(cluster.marker,
                              finalGene,by='gene')
write.csv(clusterCell,file=sprintf("%s/clusterCell.csv",FigOut))

DefaultAssay(mergData)<-'SCT'
mergData<-FindSubCluster(mergData,
                         cluster='3',
                         graph.name="SCT_snn",
                         subcluster.name='Sub3',
                         resolution=0.1)

DimPlot(mergData, reduction = "umap", group.by = c("Sub3"),label=T)

Idents(mergData)<-'Sub3'
mergData<-FindSubCluster(mergData,
                         cluster='0',
                         graph.name="SCT_snn",
                         subcluster.name='Sub0',
                         resolution=0.3)
DimPlot(mergData, reduction = "umap", group.by = c("Sub0"),label=T)
DimPlot(subset(mergData,Sub3=="0"), 
        reduction = "umap", 
        group.by = c("Sub0"),
        label=T)

Idents(mergData)<-'Sub0'
VlnPlot(mergData,
        features=c("PI16",'FGFR4'),
        raster=FALSE,
        pt.size=0) #

FeaturePlot(subset(mergData,Sub3=="0"),
            features=c('FGFR4','PI16','ACTA2','ACTG2'),
            col=c("grey","red"),
            raster=T)

DefaultAssay(mergData)<-'SCT'
Idents(mergData)<-'Sub0'
mergData<-FindSubCluster(mergData,
                         cluster='10',
                         graph.name="SCT_snn",
                         subcluster.name='Sub10',
                         resolution=0.1)

Idents(mergData)<-'Sub10'
DimPlot(mergData, reduction = "umap", group.by = c("Sub10"),label=T)
FeaturePlot(subset(mergData,Sub0=="10"),features=c('CD3D','CD8A','KLRD1'),col=c("grey","red"),raster=T)
FeaturePlot(mergData,features=c('KLRD1'),col=c("grey","red"),raster=T)
VlnPlot(mergData,features=c("KLRD1",'CD3D'),raster=FALSE,pt.size=0) #

DefaultAssay(mergData)<-'SCT'
Idents(mergData)<-'Sub10'
mergData<-FindSubCluster(mergData,
                         cluster='23',
                         graph.name="SCT_snn",
                         subcluster.name='Sub23',
                         resolution=0.1)

Idents(mergData)<-'Sub23'
DimPlot(mergData, reduction = "umap", group.by = c("Sub23"),label=T)
FeaturePlot(subset(mergData,Sub0=="23"),features=c('CD3D','CD8A','KLRD1'),col=c("grey","red"),raster=T)
FeaturePlot(mergData,features=c('KLRD1'),col=c("grey","red"),raster=T)
VlnPlot(mergData,features=c("KLRD1",'CD3D'),raster=FALSE,pt.size=0) #


DefaultAssay(mergData)<-'SCT'
Idents(mergData)<-'Sub23'
mergData<-FindSubCluster(mergData,
                         cluster='24',
                         graph.name="SCT_snn",
                         subcluster.name='Sub24',
                         resolution=0.1)
Idents(mergData)<-'Sub24'
DimPlot(mergData, reduction = "umap", group.by = c("Sub24"),label=T)
FeaturePlot(subset(mergData,Sub0=="24"),features=c('ACTA2','ACTG2','MYH11','THBS1'),col=c("grey","red"),raster=T)
VlnPlot(mergData,features=c('ACTA2','MYH11','THBS1'),raster=FALSE,pt.size=0) #

############################################################################################
############################################################################################
DefaultAssay(mergData)<-'Xenium';Idents(mergData)<-'Sub3'
DimPlot(mergData, reduction = "umap", group.by = c("ident"),label=T)
DimPlot(mergData, reduction = "umap", group.by = c("ident", "orig.ident"),label=T)
DimPlot(mergData, reduction = "umap", group.by = c("ident"),label=T)
DimPlot(mergData, reduction = "umap", group.by = c("ident"),label=T,split.by="Disease")

FeaturePlot(mergData,features=c("KRT17",'KRT5'),col=c("grey","red"),raster=T) #
VlnPlot(mergData,features=c("KRT17",'KRT5'))
FeaturePlot(mergData,features=c("AGER"),col=c("grey","red"),raster=T) #
FeaturePlot(mergData,features=c("HHIP","AGER"),col=c("grey","red"),raster=T) #
FeaturePlot(mergData,features=c("FOXJ1","CTSE"),col=c("grey","red"),raster=T) #
FeaturePlot(mergData,features=c("VWF"),col=c("grey","red"),raster=T) #
FeaturePlot(mergData,features=c("SPINK1"),col=c("grey","red"),raster=F) #
VlnPlot(mergData,features=c("SPINK1"),raster=FALSE) #
VlnPlot(mergData,features=c("KRT17","KRT5"),raster=FALSE,pt.size=0) #
VlnPlot(mergData,features=c("FOXJ1",'CTSE'),raster=FALSE,pt.size=0) #


DefaultAssay(mergData)<-'Xenium';Idents(mergData)<-'Sub0'
VlnPlot(mergData,features=c("PI16",'FGFR4'),raster=FALSE,pt.size=0) #

###############
DefaultAssay(mergData)<-'Xenium'
hiCell<-row.names(mergData@meta.data[mergData@meta.data$Sub10=='21',])
DimPlot(mergData, label=T, group.by="Sub10", raster=F,
        cells.highlight=hiCell,
        cols.highlight = c("red"), cols= "grey")+NoLegend()


DefaultAssay(mergData)<-'Xenium';Idents(mergData)<-"Sub23"
DimPlot(mergData,group.by='Sub23',label=T)
testMarker230<-FindMarkers(mergData,ident.1="23_0");head(testMarker230,20)
testMarker231<-FindMarkers(mergData,ident.1="23_1");head(testMarker231,20)
testMarker232<-FindMarkers(mergData,ident.1="23_2");head(testMarker232,20)
testMarker233<-FindMarkers(mergData,ident.1="23_3");head(testMarker233,20)

testMarker28<-FindMarkers(mergData,ident.1="28");head(testMarker28,20)
testMarker21<-FindMarkers(mergData,ident.1="21",only.pos = T);head(testMarker21,20)


DefaultAssay(mergData)<-'Xenium';Idents(mergData)<-"Sub24"
DimPlot(mergData,group.by='Sub24',label=T)
testMarker240<-FindMarkers(mergData,ident.1="24_0",only.pos = T);head(testMarker240,20)
testMarker241<-FindMarkers(mergData,ident.1="24_1",only.pos = T);head(testMarker241,20)
testMarker242<-FindMarkers(mergData,ident.1="24_2",only.pos = T);head(testMarker242,20)

FeaturePlot(mergData,features=c('ACTA2','ACTG2'),col=c("grey","red"),raster=T)
VlnPlot(mergData,features=c('ACTA2','ACTG2'),pt.size=0)

hiCell<-row.names(mergData@meta.data[mergData@meta.data$Sub10=='28',])
DimPlot(mergData, label=T, group.by="Sub10", raster=F,
        cells.highlight=hiCell,
        cols.highlight = c("red"), cols= "grey")+NoLegend()

hiCell<-row.names(mergData@meta.data[mergData@meta.data$Sub10=='24',])
DimPlot(mergData, label=T, group.by="Sub10", raster=F,
        cells.highlight=hiCell,
        cols.highlight = c("red"), cols= "grey")+NoLegend()


######for cluster 21 #####
FeaturePlot(mergData,features=c('EPCAM'),col=c("grey","red"),raster=T)
FeaturePlot(mergData,features=c('VWF'),col=c("grey","red"),raster=T)
FeaturePlot(mergData,features=c('CSF3R'),col=c("grey","red"),raster=T)


mergData$CellType<-dplyr::case_when(
  mergData$Sub24 %in% c("27")~"KRT17+KRT5-",
  mergData$Sub24 %in% c("13")~"Basal",
  mergData$Sub24 %in% c("3_1","3_2","3_3",'18')~"Ciliated",
  mergData$Sub24 %in% c("3_0")~"Club",
  mergData$Sub24 %in% c("2")~"AT2",
  mergData$Sub24 %in% c('12')~"AT1",
  
  ##FeaturePlot(mergData,features=c('MMRN1'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('22')~"Lym.En",
  ##FeaturePlot(mergData,features=c('GJA5'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('20')~"Artery.En",
  ##FeaturePlot(mergData,features=c('SHANK3','TMEM100','FCN3'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('1')~"Capillary.En",
  ##FeaturePlot(mergData,features=c('VWF','PLVAP'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('4')~"En",
  
  ##FeaturePlot(mergData,features=c('WT1','UPK3B'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('29')~"Mesothelial",
  ##FeaturePlot(mergData,features=c('FGFR4','PI16'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('0_2','0_3','0_4')~"AlvF",
  mergData$Sub24 %in% c('0_5')~"AlvF",
  ##DimPlot(subset(mergData,Sub3=="0"), reduction = "umap", group.by = c("Sub0"),label=T)
  ##FeaturePlot(subset(mergData,Sub3=="0"),features=c('ACTA2','ACTG2','MYH11','ASPN'),col=c("grey","red"),raster=T)
  ## Myo should be ACTA2+,ASPN+,ACTG2-,MYH11-
  mergData$Sub24 %in% c('0_1','24_1')~"MyoF",
  mergData$Sub24 %in% c('0_0','0_7','0_6','24_0','24_2')~"Fibro",
  ##FeaturePlot(mergData,features=c('ACTA2','ACTG2','MYH11'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('5')~"SMC",
  ##FeaturePlot(mergData,features=c('PDGFRB','NID1',"COX4I2"),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('8')~"Pericyte",
  
  
  ##FeaturePlot(mergData,features=c('MKI67','TOP2A'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('26')~"Prolif.Cells",
  ##FeaturePlot(mergData,features=c('MZB1','MS4A1','MS4A2'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('19')~"B",
  mergData$Sub24 %in% c('15')~"Plasma",
  mergData$Sub24 %in% c('16')~"Mast",
  
  #FeaturePlot(mergData,features=c('MACRO','SPP1'),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('14')~"Macro",
  mergData$Sub24 %in% c('17')~"MoM",
  #FeaturePlot(mergData,features=c("F13A1","FOLR2","IGSF21"),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('7')~"IM",
  
  #FeaturePlot(mergData,features=c('CD3E','CD40LG',"CD8A","KLRD1"),col=c("grey","red"),raster=T)
  #FeaturePlot(mergData,features=c("HAVCR2"),col=c("grey","red"),raster=T) 
  #FeaturePlot(mergData,features=c("GNLY"),col=c("grey","red"),raster=T) 
  mergData$Sub24 %in% c('6')~"CD4T",
  mergData$Sub24 %in% c('10_0','10_2')~"CD8T",
  mergData$Sub24 %in% c('10_1')~"NK",
  
  #FeaturePlot(mergData,features=c("FCGR3B","CSF3R"),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('11')~"Neutrophils",
  
  ###pDC=c("SCT","PTPRC"); DC1=c("LAMP3","CLEC9A");DC2=c("CD1E");migraDC=c("LAD1","CCL19")
  #FeaturePlot(mergData,features=c("CLEC9A",'CD1E',"LAD1","CCL19"),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('25')~"DC",
  
  #FeaturePlot(mergData,features=c("LAD1","CCL19"),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('23_0','23_1','23_2','23_3','28')~"Chemokine.Imm", 
  
  #FeaturePlot(mergData,features=c("S100A12","FCGR3A","MTSS1","CSF3R"),col=c("grey","red"),raster=T)
  #FeaturePlot(mergData,features=c("MTSS1"),col=c("grey","red"),raster=T)
  #FeaturePlot(mergData,features=c("S100A12"),col=c("grey","red"),raster=T)
  mergData$Sub24 %in% c('9')~"Monocyte",
  
  mergData$Sub24 %in% c('12')~"unknown",
  .default='unknown'
);table(mergData$Sub24,mergData$CellType)
table(mergData$Sub24,mergData$CellType)[,'unknown']
VlnPlot(mergData,features=c('nCount_Xenium', 'nFeature_Xenium'),pt.size=0) ##unknow has few cell/feature counts
FeaturePlot(mergData,features=c('CP', 'PI3','S100A7','SERPINA3'),col=c("grey","red"),raster=T)

########################################################################################################################
DefaultAssay(mergData)<-'SCT'
Idents(mergData)<-'Sub24'
mergData<-FindSubCluster(mergData,
                         cluster='19',
                         graph.name="SCT_snn", 
                         subcluster.name='Sub25',
                         resolution=0.3)
Idents(mergData)<-'Sub25'
DimPlot(mergData, reduction = "umap", group.by = c("Sub25"),label=T)
DimPlot(subset(mergData,Sub0 %in% c("19")), reduction = "umap", group.by = c("Sub25"),label=T)
FeaturePlot(subset(mergData,Sub0 %in% c("19")),features=c('MS4A1','CD8A','CD40LG'),col=c("grey","red"),raster=F)
VlnPlot(subset(mergData,Sub0 %in% c("19")),features=c('MS4A1','CD8A','CD40LG','CD3E'),raster=FALSE,pt.size=0) #
DimPlot(subset(mergData,Sub0 %in% c("19")), reduction = "umap", group.by = c("Sub25"),label=T)
FeaturePlot(subset(mergData,Sub0 %in% c("19")),features=c('MS4A1','CD8A','CD40LG','CD3E'),col=c("grey","red"),raster=F)

FeaturePlot(subset(mergData,Sub0 %in% c("19")),
            features=c('MS4A1'),
            col=c("grey","red"),
            raster=F)
FeaturePlot(subset(mergData,Sub0 %in% c("19")),
            features=c('CD79A'),
            col=c("grey","red"),
            raster=F)
FeaturePlot(subset(mergData,Sub0 %in% c("19")),
            features=c('CD8A'),
            col=c("grey","red"),
            raster=F)
FeaturePlot(subset(mergData,Sub0 %in% c("19")),
            features=c('CD40LG'),
            col=c("grey","red"),
            raster=F)
FeaturePlot(subset(mergData,Sub0 %in% c("19")),
            features=c('CD3E'),
            col=c("grey","red"),
            raster=F)

mergData$CellType<-ifelse(mergData$Sub25 %in% c("19_4"),'CD4T',mergData$CellType)
mergData$CellType<-ifelse(mergData$Sub25 %in% c("19_6"),'unknown',mergData$CellType)


saveRDS(mergData,file=sprintf("%s/IPF_mergData_annotated.RDS",saveFolder))

