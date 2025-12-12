######################
.libPaths(c("/projects/p32215/tools/forR4.4Seu5.1.0",'/home/yyw9094/R/forR4.4Seu5.1.0',"/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4"))
######################

library("Seurat")#
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
source("SpatialSourceCode.R") ##
##########################################################################################################################################
###################################                                                                 ######################################
##########################################################################################################################################

IPF_Nich<-readRDS(file=sprintf("%s/IPF_withNichAllData.RDS",saveFolder))# ###6~15 niches
table(IPF_Nich$niches12,IPF_Nich$orig.ident)
table(IPF_Nich$CellType,IPF_Nich$niches12)

#########
CTLevel<-c('AT2','AT1','Ciliated','Club','Basal','KRT17+KRT5-',
           'En',"Capillary.En","Artery.En","Lym.En",
           "Fibro","AlvF","MyoF","SMC","Pericyte","Mesothelial",
           "B","Plasma","CD4T","CD8T","NK",
           "Monocyte","DC","Macro","MoM","IM","Neutrophils","Chemokine.Imm","Prolif.Cells","Mast"
           );length(CTLevel)
IPF_Nich$CellType<-factor(IPF_Nich$CellType,levels=(CTLevel))
Idents(IPF_Nich)<-"CellType"

#########
library(RColorBrewer)
require(gridExtra)
colors_df <- colorRampPalette(brewer.pal(9, "Set1"))(length(unique(IPF_Nich$CellType)));length(colors_df)
colors_df<-c("#E41A1C", "#B43547", "#845172", "#546C9D", "#3983AC", "#3F908E", "#FFB315", "#FFD723", "#FFFA31", "#E9D630",
             "#85658D", "#9B4F9D", "#B75D70", "#D46A43", "#F07816", "#FF9007",
             'blue','cyan','green1','green4',####. B-CD8T   "#FFB315", "#FFD723", "#FFFA31", "#E9D630",
             "#D0A72D", "#B8782A", "#AB5832", "#C1645C",
             "#D87085", "#EE7CAF", "#E685B8", "#CC8BAD", "#B292A3", "#999999");length(colors_df)

############################################################################################################
############################################################################################################
############################################################################################################
library(ggplot2)
library(ggalluvial)
library(dplyr)

colors_Nich <- colorRampPalette(brewer.pal(12, "Set1"))(12);length(colors_Nich)
colors_Nich<-c("blue", "#419486", "#5A9D5A", "#91569A","#E41A1C", "orange", "coral", "#F6EF32", "#B6742A",'cyan', "#DD87B4", "#999999")
nichSel<-paste("Niche",1:15,sep="");nichSel

head(IPF_Nich)
plt_Prop_Allu<-function(datIn=IPF_Nich@meta.data,
                        YBar='niches12',
                        xCat='CellType',
                        rowNPaste='Niche',
                        xLevel=CTLevel,
                        aluWidth=0.5,
                        colSet='Set1',
                        colSel=NULL,
                        nichSel=nichSel,
                        colLen=9
                        ){
  tabIn<-table(datIn[,YBar],datIn[,xCat])
  Cell_prop <-as.data.frame.matrix(sweep(tabIn, 2, colSums(tabIn), FUN = "/"))
  if(!is.null(rowNPaste)){
    row.names(Cell_prop)<-paste(rowNPaste,row.names(Cell_prop),sep="")
  }
  Cell_prop$RowName<-row.names(Cell_prop)
  
  longDf<-tidyr::pivot_longer(Cell_prop,  cols = -RowName,  names_to = "CT", values_to = "Prop")
  pltData<-longDf
  
  NicheLevel<-row.names(Cell_prop)
  if (!is.null(colSel)){
    colIn<-colSel
  }else{
    colIn<-colorRampPalette(brewer.pal(colLen, colSet))(length(NicheLevel))
  }
  
  colSch<-data.frame(NicheLevel=NicheLevel, col=colIn)
  if(!is.null(nichSel)){
    tempNichL<-nichSel[nichSel %in% colSch$NicheLevel]
    colSch$NicheLevel<-factor(colSch$NicheLevel,levels=tempNichL)
    colSch<-colSch[order(colSch$NicheLevel),]
  }
  pltData$xCat<-factor(pltData$RowName,levels=colSch$NicheLevel)
  pltData$CT<-factor(pltData$CT,levels=xLevel)
  
  (aluPlot<-ggplot(pltData, aes(x = CT, 
                                stratum = xCat, 
                                alluvium = xCat,
                                y = Prop)) +
      geom_alluvium(aes(fill = xCat), width = aluWidth) +
      geom_stratum(width = aluWidth, aes(fill = xCat)) +
      scale_fill_manual(values = colSch$col)+
      theme_minimal() +
      labs(title = "", x = "", y = "Proportion") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)))

  return(aluPlot)
}

allNich<-c("niches12")
lapply(1:length(allNich),function(iX){#
  inNich<-allNich[iX]
  colors_NichIn=colors_Nich
  if (!iX==7){colors_NichIn=NULL}
  pl<-plt_Prop_Allu(datIn=IPF_Nich@meta.data,
                    YBar=inNich,
                    xCat='CellType',
                    rowNPaste='Niche',
                    xLevel=CTLevel,
                    aluWidth=0.6,
                    colSet='Paired',
                    colLen=12,
                    colSel=colors_NichIn,
                    nichSel=nichSel)
  
  pdf(sprintf("%s/Compostion_%s.pdf",FigOut, inNich),width=7.87,height=7.87/1.5, family="ArialMT")
  print(pl)
  dev.off()  
})


############################################################################################################
############################################################################################################
############################################################################################################
IPF_Nich$NichF<-IPF_Nich$niches12;table(IPF_Nich$NichF)
IPF_Nich$NichF<-factor(IPF_Nich$NichF,levels=1:12)

IPF_Nich$FNich<-case_when(
  IPF_Nich$NichF %in% c('1')~'TLS',
  IPF_Nich$NichF %in% c('2')~'MoM',
  IPF_Nich$NichF %in% c('3')~'AM',
  IPF_Nich$NichF %in% c('4')~'SMC',
  IPF_Nich$NichF %in% c('5')~'B_Zone',
  IPF_Nich$NichF %in% c('6')~'Alveolar',
  IPF_Nich$NichF %in% c('7')~'Airway',
  IPF_Nich$NichF %in% c('8')~'Plasma',
  IPF_Nich$NichF %in% c('9')~'Mixed1',
  IPF_Nich$NichF %in% c('10')~'FF',
  IPF_Nich$NichF %in% c('11')~'Mixed2',
  IPF_Nich$NichF %in% c('12')~'Artery',
  .default='unknown'
  );table(IPF_Nich$FNich,IPF_Nich$NichF)
IPF_Nich$FNich<-factor(IPF_Nich$FNich,levels=c("TLS",'MoM','AM','SMC','TLS_GC','Alveolar','Airway','Plasma','Mixed1','FF','Mixed2','Artery'))


plot(1:12,1:12,col=colors_Nich)
names(colors_Nich)<-levels(IPF_Nich$FNich)

(uniID<-unique(IPF_Nich$orig.ident))
lapply(1:length(uniID),function(iD){
  imp<-ImageDimPlot(IPF_Nich, 
                    fov = sprintf('fov_%s',uniID[iD]),
                    group.by = "FNich", 
                    size = 0.5, 
                    dark.background = FALSE) + ggtitle("Niches") +
    scale_fill_manual(values =colors_Nich)
  print(uniID[iD])
  
  sapply(1:2,function(ix){
    if(ix==1){pdf(sprintf("%s/Nich_%s.pdf",FigOut,uniID[iD]),width=7.87,height=7.87)}
    if(ix==2){tiff (sprintf("%s/Nich_%s.tiff",FigOut,uniID[iD]), res=600, width=7.87,height=7.87,units="in",compression="lzw")}
    print(imp)
    dev.off()})
})


names(colors_df)<-levels(IPF_Nich$CellType)
lapply(1:length(uniID),function(iD){
  imp<-ImageDimPlot(IPF_Nich, 
                    fov = sprintf('fov_%s',uniID[iD]),
                    group.by = "CellType", 
                    size = 0.5, 
                    dark.background = FALSE) + ggtitle("CellType") +
    scale_fill_manual(values =colors_df)
  print(uniID[iD])
  
  sapply(1:2,function(ix){
    if(ix==1){pdf(sprintf("%s/AllCT_%s.pdf",FigOut,uniID[iD]),width=7.87*1.2,height=7.87)}
    if(ix==2){tiff (sprintf("%s/AllCT_%s.tiff",FigOut,uniID[iD]), res=600, width=7.87,height=7.87,units="in",compression="lzw")}
    print(imp)
    dev.off()})
})


#####################
IPF_Nich$ptAll<-ifelse(IPF_Nich$orig.ident %in% c("IPF3L",'IPF3R'), 'IPF3', IPF_Nich$orig.ident)
(ptNich<-table(IPF_Nich$ptAll,IPF_Nich$FNich))
prop.table(ptNich, margin=1)
propNich<-as.data.frame(prop.table(ptNich, margin=1));head(propNich)
colnames(propNich)<-c("PT",'Nich','Prop')

(pltniCop<-ggplot(propNich, aes(x = PT, 
                    stratum = Nich, 
                    alluvium = Nich,
                    y = Prop)) +
  geom_alluvium(aes(fill = Nich), width = 0.6) +
  geom_stratum(width = 0.6, aes(fill = Nich)) +
  scale_fill_manual(values = colors_Nich)+
  theme_minimal() +
  labs(title = "", x = "", y = "Proportion") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))

pdf(sprintf("%s/Nich_ComEachPT.pdf",FigOut),width=7.87*1.5,height=7.87/1.2)
print(pltniCop)
dev.off()
