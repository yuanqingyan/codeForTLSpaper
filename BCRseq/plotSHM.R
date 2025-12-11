######################
.libPaths(c("/projects/p32215/tools/forSlurmR440","/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4","/projects/p32215/tools/forR4.4Seu5.1.0"))
######################

###########################################################################
######### Note: the table_SHM.csv files were exported from platforma    ###
######### We used platforma for the analysis, including clotyping et al. ##
###########################################################################

SaveFolder<-"./ManuscriptData/ILD_B/BCR"
list.files(SaveFolder)
list.files(sprintf("%s/SHM",SaveFolder))

FigOut<-"/projects/b1042/Yuanqing/RNAseq/BCR"
setwd(FigOut)

#############
ctrlID<-c('Donor_9','Donor_4','Donor_2','Donor_7', 'Donor_5','Donor_6');length(ctrlID)

Light<-read.delim(sprintf("%s/SHM/L_table_SHM.csv",SaveFolder),header=T,sep=',',row.names=1)
sum(ctrlID %in% Light$Sample)

colnames(Light)[c(1:9,21)]<-c('TreeID','Sample','Nnodes','Nclone','Mmole','TotalReadCount','Chain','V','J','CDR3MRA');head(Light)
Light$indx<-1:nrow(Light)
Light[is.na(Light$Nclone),]
Light$disease<-ifelse(Light$Sample %in% ctrlID,'Ctrl','PF')
Light$disease<-factor(Light$disease,levels=c("PF",'Ctrl'))

library(ggplot2)
p_H<-ggplot(Light, aes(x = indx, y = Nclone)) +
  geom_point(aes(color = disease)) +
  facet_wrap(~ disease, ncol = 2) +
  labs(title = "Nclone", x = "indx",y = "Nclone") +
  ggtitle('Number of clone per tree')+
        theme_bw()+
        theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.5,size=12,colour = "black"),
              axis.title=element_text(size=12,colour = "black"),
              axis.text.y = element_text(size=12,colour = "black"),
              panel.border = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              legend.position = "none")

pdf(sprintf("%s/Ncl_perTree_L.pdf",FigOut),width=7.87,height=7.87/2, family="ArialMT")
print(p_H)
dev.off()




#####################################################################
Heavy<-read.delim(sprintf("%s/SHM/H_table_SHM.csv",SaveFolder),header=F,sep=',',row.names=1)
colnames(Heavy)<-c('TreeID','Sample','Nnodes','Nclone','Mmole','TotalReadCount','Chain','V','J','CDR3MRA')
Heavy$indx<-1:nrow(Heavy)
Heavy[is.na(Heavy$Nclone),]
Heavy$disease<-ifelse(Heavy$Sample %in% ctrlID,'Ctrl','PF')
table(Heavy$disease) 
Heavy$disease<-factor(Heavy$disease,levels=c("PF",'Ctrl'))

library(ggplot2)
p_H<-ggplot(Heavy, aes(x = indx, y = Nclone)) +
  geom_point(aes(color = disease)) +
  facet_wrap(~ disease, ncol = 2) +
  labs(title = "Nclone", x = "indx",y = "Nclone") +
  ggtitle('Number of clone per tree')+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.5,size=12,colour = "black"),
        axis.title=element_text(size=12,colour = "black"),
        axis.text.y = element_text(size=12,colour = "black"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none")

pdf(sprintf("%s/Ncl_perTree.pdf",FigOut),width=7.87,height=7.87/2, family="ArialMT")
print(p_H)
dev.off()

#####################################################################
Heavy$Vgene<-sapply(strsplit(Heavy$V,split="\\-"),function(x) x[1])
Heavy_sub<-Heavy[Heavy$Nclone>=25,]
table(Heavy_sub$disease)


pdf(sprintf("%s/barplot_25Cut_Vgene.pdf",FigOut),width=7.87/1.8,height=7.87/2, family="ArialMT")
print(barplot(table(Heavy_sub$Vgene),col=c('grey','lightblue','tomato','pink','lightgreen'),las=2))
dev.off()

