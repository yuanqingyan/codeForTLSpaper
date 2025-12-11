######################
.libPaths(c("/projects/p32215/tools/forSlurmR440","/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4","/projects/p32215/tools/forR4.4Seu5.1.0"))
######################



generateSubCl<-function(obj=seuST,proN="AdvF",nFeature=3000,npcs=35,res=0.9,kparam=50){
  temp_CT<-obj
  rawCount_temp<-temp_CT@assays$RNA@counts
  har_temp<-CreateSeuratObject(counts = rawCount_temp, project =proN, min.cells = 3) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = nFeature) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(npcs = npcs, verbose = FALSE)
  har_temp<-AddMetaData(har_temp,metadata=temp_CT@meta.data[,c(4:ncol(temp_CT@meta.data))])
  har_temp<-har_temp %>% RunHarmony("sampleID", plot_convergence = FALSE)
  har_temp <- har_temp %>%
    RunUMAP(reduction = "harmony", dims = 1:npcs) %>%
    FindNeighbors(reduction = "harmony", dims = 1:npcs, k.param = kparam) %>%
    FindClusters(resolution = res) %>%
    identity()
  return(har_temp)
}
########################################

library(Seurat)
library(harmony)
library(cowplot)
library(dplyr)
library(ggplot2)
FigOut<-"/projects/b1042/Yuanqing/SC/ILD/output"


##########################################################################
SunB<-readRDS(file="./scRNAseq_B_Plasma_SunDataset.RDS")
##########################################################################

DimPlot(SunB,group.by='Bsub',label=T,repel=T)
DimPlot(SunB,group.by='Bsub',label=T,repel=T,split.by='whichDisease',ncol=2)+NoLegend()

gcbGene<-c("SERPINA9","MEF2B","RGS13")
navGene<-c("TCL1A",'CD72','IGHD','IGHM')
mebGene<-c("CD27",'CD69')
FeaturePlot(SunB,features=gcbGene,col=c('grey','red'))
FeaturePlot(SunB,features=navGene,col=c('grey','red'))

FeaturePlot(SunB,features=c("PAICS",'BZW2','ABCE1'),col=c('grey','red'))
FeaturePlot(SunB,features=c("CD27",'CD69'),col=c('grey','red'))



SunB2<-generateSubCl(obj=SunB,
                      proN='B2',
                      nFeature=3000,
                      npcs=35,
                      res=1,
                      kparam=30)
head(SunB2)
DefaultAssay(SunB2)<-"RNA";Idents(SunB2)<-'RNA_snn_res.1'
DimPlot(SunB2,group.by='RNA_snn_res.1',label=T)+NoLegend()
DimPlot(SunB2,group.by='Bsub',label=T,repel=T)+NoLegend()

FeaturePlot(SunB2,features=gcbGene,col=c('grey','red'))
FeaturePlot(SunB2,features=navGene,col=c('grey','red'))
FeaturePlot(SunB2,features=mebGene,col=c('grey','red'))
FeaturePlot(SunB2,features=c("PAICS",'BZW2','ABCE1'),col=c('grey','red'))
FeaturePlot(SunB2,features=c("MS4A1",'IGHA1','IGKC'),col=c('grey','red'))
FeaturePlot(SunB2,features=c("MS4A1"),col=c('grey','red'))
FeaturePlot(SunB2,features=c("MZB1"),col=c('grey','red'))
FeaturePlot(SunB2,features=c('MKI67'),col=c('grey','red'))

##########
SunB2<-FindSubCluster(
  SunB2,
  cluster='7',
  graph.name='RNA_snn',
  subcluster.name = "sub7",
  resolution = 0.3,
  algorithm = 1
)
DimPlot(SunB2, group.by='sub7',label=T, repel=T,label.size=4)+NoLegend()

Idents(SunB2)<-'sub7'
SunB2<-FindSubCluster(
  SunB2,
  cluster='7_1',
  graph.name='RNA_snn',
  subcluster.name = "sub71",
  resolution = 0.9,
  algorithm = 1
)
DimPlot(SunB2, group.by='sub7',label=T, repel=T,label.size=4)+NoLegend()
Idents(SunB2)<-'sub71'
VlnPlot(SunB2,features=c('MKI67','MS4A1','IGKC'),pt.size=0)

FeaturePlot(SunB2,features=c('MS4A1'),col=c('grey','red'))
FeaturePlot(SunB2,features=c('MZB1'),col=c('grey','red'))
FeaturePlot(SunB2,features=c('IGKC'),col=c('grey','red'))


DefaultAssay(SunB2)<-'RNA';Idents(SunB2)<-"sub7"
sunMark<-FindAllMarkers(SunB2,min.pct=0.1,only.pos=T)
sunMark<-sunMark[order(sunMark$cluster,-sunMark$avg_log2FC),]#;head(sunMark4,20)
sunMark %>% group_by(cluster) %>% top_n(10, avg_log2FC) %>% as.data.frame()

SunB2$FineCT<-dplyr::case_when(
  SunB2$sub7 %in% c('7_1')~'Prolif.B',
  SunB2$sub7 %in% c('7_0')~'GCB',
  SunB2$sub7 %in% c('5')~'Naive.B',
  SunB2$sub7 %in% c('0','3')~'Memory.Act.B',
  .default='Plasma'
);DimPlot(SunB2,group.by='FineCT',label=T,repel=T)+NoLegend()

DimPlot(SunB2,group.by='FineCT',split.by="whichDisease")


saveRDS(SunB2,file="./SunDataSet_FineCT.RDS")
############
SunB2<-readRDS(file="./SunDataSet_FineCT.RDS")
SunB2$FineCT<-factor(SunB2$FineCT,levels=c("Naive.B",
                                               "Memory.Act.B",
                                               "GCB",
                                               "Prolif.B",
                                               'Plasma'))
colPale<-c("#00FFFF", "#800080", "#FF0000", "#0000FF", "#EE82EE")



pdf(sprintf("%s/Sun.Dimplot.pdf",FigOut),width=7.87*1.2,height=7.87, family="ArialMT")
print(DimPlot(SunB2,group.by="FineCT",label=F,repel=T,cols=colPale) & NoAxes())
dev.off()



library(ggplot2)
pltMarkerGene<-c("MS4A1",
                 "TCL1A","CD72","IGHD","IGHM",
                 "CD69", "CD27",'TNFSF9','CR2',
                 "SERPINA9", "MEF2B", "RGS13",
                 "MKI67",'TOP2A',
                 "MZB1",'IGKC','IGHA1','JCHAIN')

SunB2$FineCT<-factor(SunB2$FineCT,levels=rev(c("Naive.B",
                                                   "Memory.Act.B",
                                                   "GCB",
                                                   "Prolif.B",
                                                   'Plasma')))

Idents(SunB2)<-"FineCT"
dotT<-DotPlot(object =SunB2, features = pltMarkerGene,cols = c("lightgrey", "red"))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dotT

pdf(sprintf("%s/Sun.Dotplot.pdf",FigOut),width=7.87*1.2,height=7.87/1.2, family="ArialMT")
print(dotT)
dev.off()


SunB2$FineCT<-factor(SunB2$FineCT,levels=c("Naive.B",
                                           "Memory.Act.B",
                                           "GCB",
                                           "Prolif.B",
                                           'Plasma'))

pdf(sprintf("%s/Sun.Dimplot_split.pdf",FigOut),width=7.87*1.2*1.5,height=7.87/2, family="ArialMT")
print(DimPlot(SunB2,
              group.by="FineCT",
              label=F,
              repel=T,
              cols=colPale,
              split.by='whichDisease') & NoAxes())
dev.off()



