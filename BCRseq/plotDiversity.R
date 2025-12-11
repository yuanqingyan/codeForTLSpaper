######################
.libPaths(c("/projects/p32215/tools/forSlurmR440","/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4","/projects/p32215/tools/forR4.4Seu5.1.0"))
######################


###########################################################################
######### Note: the table_SHM.csv files were exported from platforma    ###
######### We used platforma for the analysis, including clotyping et al. ##
###########################################################################


SaveFolder<-"./ManuscriptData/ILD_B/BCR"
list.files(sprintf("%s/Diversity",SaveFolder),pattern=".csv")
FigOut<-"/projects/b1042/Yuanqing/RNAseq/BCR"

Heavy<-read.csv(sprintf("%s/Diversity/RepertoireDiversity_Heavt.csv",SaveFolder),row.names=2)
Light<-read.csv(sprintf("%s/Diversity/RepertoireDiversity_Light.csv",SaveFolder),row.names=2)
unique(row.names(Heavy))
contrlID<-c("Donor_9","Donor_4","Donor_2","Donor_7","Donor_5","Donor_6");length(contrlID)

##Observed.Diversity → Chao1.Estimate → Efron.Thisted.Estimate → Shannon.Wiener.Diversity / Shannon.Wiener.Index → Normalized.Shannon.Wiener.Index → Inverse.Simpson.Index → Gini.Index → D50.Diversity
RichnessToEvenness<-c("Observed.Diversity", "Chao1.Estimate", "Efron.Thisted.Estimate", "Shannon.Wiener.Diversity", "Shannon.Wiener.Index",
                      "Normalized.Shannon.Wiener.Index", "Inverse.Simpson.Index", "Gini.Index", "D50.Diversity")
Heavy<-Heavy[,c("X.",RichnessToEvenness)]
Light<-Light[,c("X.",RichnessToEvenness)]

library(coin)
allDiv<-list(Heavy=Heavy, Light=Light)
allStat_fromCoin<-lapply(1:length(allDiv), function(x){
  temp<-allDiv[[x]]
  colnames(temp)[1]<-'grp'
  temp$grp<-ifelse(row.names(temp) %in% contrlID, 'Control','PF')
  allP<-as.data.frame(do.call('rbind',lapply(2:ncol(temp),function(y){
    temp2<-temp[,c(1,y)]
    colnames(temp2)[2]<-'Diversity'
    temp2$grp<-factor(temp2$grp, levels=c("Control",'PF'))
    p<-pvalue(wilcox_test(Diversity~grp, data=temp2))#$p.value
    return(p)
  })));allP$Index<-colnames(temp)[2:ncol(temp)]
  allP$fdr<-p.adjust(allP$V1,'BH')
  return(allP)
});allStat_fromCoin


####################################################
library(ggplot2)
allplt<-lapply(1:length(allDiv), function(x){
  temp<-allDiv[[x]]
  colnames(temp)[1]<-'grp'
  temp$grp<-ifelse(row.names(temp) %in% contrlID, 'Control','PF')
  subPlot<-lapply(2:ncol(temp),function(y){
    temp2<-temp[,c(1,y)]
    colnames(temp2)[2]<-'Diversity'
    temp2$grp<-factor(temp2$grp, levels=c("Control",'PF'))
    
    bp_col<-c('green','magenta')
    
    set.seed(123)
    p <- ggplot(temp2, aes(x=grp, y=Diversity,col=grp)) + 
      geom_boxplot()+
      geom_jitter(shape=1, position=position_jitter(0.3),size=3)+
      scale_color_manual(values=bp_col)+
      labs(y = "Index")+
      labs(x = "")+
      ggtitle(sprintf("%s",colnames(temp)[y]))+
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
  })
  return(subPlot)
});allplt

length(allplt[[1]])
cplt<-c(allplt[[1]],allplt[[2]])
require(gridExtra)
pdf(sprintf("%s/H_L_Diversity.pdf",FigOut),width=7.87*3,height=7.87, family="ArialMT")
print(do.call(grid.arrange, c(cplt, ncol=9,nrow=2)))
dev.off()

