######################
.libPaths(c("/projects/p32215/tools/forSlurmR440","/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4","/projects/p32215/tools/forR4.4Seu5.1.0"))
######################

###############################################################################
###Note: 'ILD_XXX' is changed to 'PF_XXX' to make consistence of the paper ####
###############################################################################


library(limma)
library(dplyr)
library(gplots)
library(msigdbr)
library(fgsea)
library(ggplot2)
library(RColorBrewer)
require(gridExtra)



saveFolder="./TMT/RDS"

norDat<-read.csv(sprintf('%s/tc-1251_proteinquant_3948_protein.csv',saveFolder),header=T,skip=1)
row.names(norDat)<-norDat$ProteinID
sum(duplicated(norDat$GeneSymbol))
sum(duplicated(norDat$ProteinID))
table(norDat$Peptides)

norDat$sd<-apply(norDat[,which(colnames(norDat)=="Ctrl_1"):which(colnames(norDat)=="ILD_6")],1,function(x) sd(x))
norDat$mean<-rowMeans(norDat[,which(colnames(norDat)=="Ctrl_1"):which(colnames(norDat)=="ILD_6")])
norDat$CV<-norDat$sd/norDat$mean
norDat<-norDat[order(-norDat$CV),]

inDat<-norDat[,which(colnames(norDat)=="Ctrl_1"):which(colnames(norDat)=="ILD_6")]

(group<-sapply(strsplit(colnames(inDat),split="\\_"),function(x) x[1]))
(design <- model.matrix(~ group))

fit <- lmFit(inDat, design, method='robust', maxit=1000)
fit <- eBayes(fit,robust=TRUE)


allStat<-topTable(fit, coef=2, adjust="BH",n = Inf)
plot(allStat$logFC,-log10(allStat$adj.P.Val))

allStat$ProteinID<-row.names(allStat)
allStat1<-dplyr::left_join(allStat,norDat[,c('ProteinID',"GeneSymbol","Description",'Peptides')],by=c("ProteinID"))
sum(allStat1$adj.P.Val<0.05)


sum(duplicated(allStat1$GeneSymbol))
allStat1<-allStat1[order(allStat1$P.Value),]
#####
dupStat<-allStat1[allStat1$GeneSymbol %in% allStat1$GeneSymbol[duplicated(allStat1$GeneSymbol)],]
dupStat<-dupStat[order(dupStat$GeneSymbol),]

######################################################
sigPro<-allStat1[allStat1$adj.P.Val<0.05 & abs(allStat1$logFC)>log2(2),];dim(sigPro)
write.csv(sigPro,file=sprintf("%s/sigProtein.csv",saveFolder))

selStat_sig<-allStat1[allStat1$ProteinID %in% sigPro$ProteinID,c("ProteinID", "GeneSymbol")]
plotDat0_sig<-inDat[row.names(inDat) %in% selStat_sig$ProteinID,]
plotDat0_sig$ProteinID<-row.names(plotDat0_sig)
plotDat1_sig<-dplyr::left_join(plotDat0_sig,selStat_sig,by="ProteinID")
row.names(plotDat1_sig)<-plotDat1_sig$ProteinID


##########################################################
#####  Gene Set Enrichment analysis  #####################
##########################################################
allStat1<-allStat1[order(allStat1$P.Value),]
stat_forGSEA<-allStat1[!duplicated(allStat1$GeneSymbol),]
ranks <- stat_forGSEA$t
names(ranks) <- stat_forGSEA$GeneSymbol
ranks<-ranks[order(-ranks)]

stats <- ranks
geneSets <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = "GO:BP")
geneSets <- geneSets[geneSets$gene_symbol %in% names(ranks),]
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
eaRes <- fgsea(pathways = m_list, stats = stats, nperm = 5e4, minSize = 10)
eaRes <-eaRes[order(eaRes$padj),]
fwrite(eaRes,file=sprintf("%s/GOBP_GSEA.csv",saveFolder))

#######################################################################
allData<-list(
  group=group,
  allStat1=allStat1,
  sigPro=sigPro,
  selStat_sig=selStat_sig,

  stats=stats,
  eaRes=eaRes,
  m_list=m_list,

  norDat=norDat,
  inDat=inDat)

saveRDS(allData,file=sprintf("%s/All_TMT_Data.RDS",saveFolder))
#######################################################################


#######################################################################
allData<-readRDS(file=sprintf("%s/All_TMT_Data.RDS",saveFolder))

group=allData$group
allStat1=allData$allStat1
sigPro=allData$sigPro
selStat_sig=allData$selStat_sig

stats=allData$stats
eaRes=allData$eaRes
m_list=allData$m_list

norDat=allData$norDat
inDat_f=allData$inDat

#######################################################################
#######     volcano plot              #################################
#######################################################################
FigOut="/projects/b1042/Yuanqing"

pdf(sprintf("%s/Volcano.pdf",FigOut),width=6,height=6)
allStat1$col<-ifelse(allStat1$adj.P.Val<0.05,'red','black')
plot(allStat1$logFC, 
     -log10(allStat1$adj.P.Val),
     xlab='LogFC',
     ylab="-log10 (FDR)",
     col=allStat1$col)
dev.off()

#######################################################################
#######  plot each individual protein #################################
#######################################################################

bp_proteinPlt<-function(proteinID="Q9Y2Y8", 
                        colo=c("light blue",'magenta')){
  
  proteinName=norDat[proteinID,'GeneSymbol']
  df.plot<-data.frame(exp=as.numeric(inDat_f[proteinID,]),group=group)
  bp_col<-colo
  
  set.seed(123)
  p <- ggplot(df.plot, aes(x=group, y=exp,col=group)) + 
    geom_boxplot()+
    geom_jitter(shape=1, position=position_jitter(0.3),size=3)+ #,colour = df.plot$jitCol
    scale_color_manual(values=bp_col)+
    labs(y = "Expression Level")+
    labs(x = "")+
    ggtitle(sprintf("%s",proteinName))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.5,size=12,colour = "black"),
          axis.title=element_text(size=12,colour = "black"),
          axis.text.y = element_text(size=12,colour = "black"),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position = "none") 
  return(p)
}
  
width=3; height=3

allProteinID<-c("Q9Y2Y8",'P01876','P01591')
allStat1[allStat1$ProteinID %in% allProteinID,]
lapply(allProteinID,function(x){
  proteinName=norDat[x,'GeneSymbol']
  pdf(sprintf("%s/Protein_BoxPlt_%s_%s.pdf", FigOut,x,proteinName), width=width,height=height)
  p<-bp_proteinPlt(proteinID=x,colo=c("blue",'magenta'))
  print(p)
  dev.off()
})


#######################################################################
#######     dot plot of GOBP          #################################
#######################################################################
bp<-eaRes[order(-eaRes$NES),]
bp<-bp[1:15,]
bp$Description<-bp$pathway
bp$qvalue<-bp$padj
bp$p.adjust<-bp$padj
bp$FoldEnrichment<-bp$NES
bp$Count<-unlist(sapply(strsplit(as.character(bp$leadingEdge),split="\\,"),function(x) {length(x)}))


eaRes_for_plot <- bp %>%
  dplyr::arrange(padj) %>%  # Sort by significance
  dplyr::mutate(log_padj = -log10(padj)) %>%  # Convert p-value for visualization
  head(15)  # Top pathways
eaRes_for_plot$shortP<-gsub("GOBP_","",eaRes_for_plot$pathway)
eaRes_for_plot$shortP<-tolower(gsub("_"," ",eaRes_for_plot$shortP))

pltBP<-ggplot(eaRes_for_plot, 
              aes(x = NES, y = reorder(shortP, NES), size = 1, color = padj)) +
  geom_point() +
  scale_color_gradient(low = "red", high = "darkred") +
  theme_bw()+
  labs(x = "NES", y = "Pathway", 
       size = "FDR", color = "FDR",
       title = "")+
  theme(axis.text.x = element_text(angle = 00,hjust=0.95,vjust=0.5,size=12,colour = "black"),
        axis.title=element_text(size=12,colour = "black"),
        axis.text.y = element_text(size=12,colour = "black"),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) 

pdf(sprintf("%s/BP_DotPlot.pdf",FigOut), width=7,height=6)
print(pltBP)
dev.off()


#######################################################################
#######     plot select pathways      #################################
#######################################################################
selPath<-c("GOBP_GLYCOSYLATION",
           "GOBP_POSITIVE_REGULATION_OF_B_CELL_ACTIVATION",
           'GOBP_IMMUNOGLOBULIN_PRODUCTION',
           'GOBP_REGULATION_OF_B_CELL_ACTIVATION',
           'GOBP_B_CELL_PROLIFERATION')

lapply(1:length(selPath),function(x){
  p1<-plotEnrichment(m_list[[selPath[x]]], stats) + labs(title=selPath[x])
  
  pdf(sprintf("%s/%s.pdf",FigOut,selPath[x]), width=6,height=4)
  print(p1)
  dev.off()
  
  ###############################################################
  (selGene<-eaRes[eaRes$pathway %in% selPath[x],]$leadingEdge[[1]])
  
  selStat<-allStat1[allStat1$GeneSymbol %in% selGene,c("ProteinID", "GeneSymbol")]
  plotDat0<-inDat[row.names(inDat) %in% selStat$ProteinID,]
  plotDat0$ProteinID<-row.names(plotDat0)
  plotDat1<-dplyr::left_join(plotDat0,selStat,by="ProteinID")
  row.names(plotDat1)<-paste(plotDat1$GeneSymbol,plotDat1$ProteinID,sep="__")
  
  plotDat<-plotDat1[,1:11]
  color.map<-rep(c("green","magenta"),c(5,6))
  mini<-min(plotDat);max<-max(plotDat)
  len_bk=100
  bk<-seq(mini,max,by=(max-mini)/len_bk)
  mycols<-c(colorRampPalette(colors = c("darkblue","blue","white"))(length(bk)),
            colorRampPalette(colors = c("white","red","darkred"))(length(bk)))
  dist.my <- function(x) as.dist(1-cor(t(x)))
  hclust.my <- function(x) hclust(x, method="complete")
  
  pdf(sprintf("%s/heat_%s.pdf",FigOut,selPath[x]),width=10,height=10)
  sidebarcolors <- color.map
  heatmap.2(as.matrix(plotDat),
            Rowv=TRUE,
            Colv=FALSE,
            cexRow=1.2,
            cexCol=1.2,
            distfun=dist.my,
            hclustfun = hclust.my,
            ColSideColors=sidebarcolors,
            dendrogram=c("row"),
            col=mycols,
            key=TRUE,
            keysize=0.6,
            symkey=TRUE,
            scale="row",
            density.info="none",
            trace="none",
            margins=c(8,14),
            main="")
  dev.off()
  
})
