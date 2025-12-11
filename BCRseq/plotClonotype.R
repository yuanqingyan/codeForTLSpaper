######################
.libPaths(c("/projects/p32215/tools/forSlurmR440","/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4","/projects/p32215/tools/forR4.4Seu5.1.0"))
######################


###########################################################################
######### Note: the table_SHM.csv files were exported from platforma    ###
######### We used platforma for the analysis, including clotyping et al. ##
###########################################################################


SaveFolder<-"./ManuscriptData/ILD_B/BCR"
list.files(SaveFolder)
list.files(sprintf("%s/Clonotype",SaveFolder))
FigOut<-"/projects/b1042/Yuanqing/RNAseq/BCR"
setwd(FigOut)

Heavy<-read.delim(sprintf("%s/Clonotype/clones.tsv",SaveFolder),header=T,sep='\t',row.names=1)
Light<-read.delim(sprintf("%s/Clonotype/clones_light.tsv",SaveFolder),header=T,sep='\t',row.names=1)

#########################################################
getTargCol<-function(targetCol=targetCol, chain=Heavy){
  tempMat<-chain[,grepl("^Fraction.of.UMIs",colnames(chain))]
  colnames(tempMat)<-sapply(strsplit(colnames(tempMat),split="...",fix=TRUE),function(x) x[2])
  
  tempMat$targetColumn<-chain[,targetCol]
  list_temp_target<-lapply(1:(ncol(tempMat)-1),function(x) {
    temp2<-tempMat[,c(x,ncol(tempMat))]
    temp2<-temp2[!is.na(temp2[,1]),]
    return(temp2)
  });names(list_temp_target)<-colnames(tempMat)[1:(ncol(tempMat)-1)]
  return(list_temp_target)
}


list_H_CDR3Length<-getTargCol(targetCol='Length.of.CDR3.aa', chain=Heavy)
list_L_CDR3Length<-getTargCol(targetCol='Length.of.CDR3.aa', chain=Light)
list_H_BestV<-getTargCol(targetCol='Best.V.gene', chain=Heavy)
list_L_BestV<-getTargCol(targetCol='Best.V.gene', chain=Light)


#############################################################################################
summary_CDR3Length_H<-lapply(list_H_CDR3Length,function(x){
  mean<-mean(x$targetColumn)
  median<-median(x$targetColumn)
  q1<-as.numeric(quantile(x$targetColumn,0.1))
  q90<-as.numeric(quantile(x$targetColumn,0.75))
  return(list(mean=mean,median=median,q1=q1,q90=q90))
})

df_H_CDR3<-as.data.frame(unlist(lapply(summary_CDR3Length_H,function(x) x$q90)))
colnames(df_H_CDR3)<-'mean'
df_H_CDR3$grp<-ifelse(row.names(df_H_CDR3) %in% c('Donor_9','Donor_4','Donor_2','Donor_7', 'Donor_5','Donor_6'),'Control','PF')
df_H_CDR3$grp<-as.factor(df_H_CDR3$grp)
boxplot(mean~grp,data=df_H_CDR3)
t.test(mean~grp,data=df_H_CDR3)
coin::wilcox_test(mean~grp,data=df_H_CDR3)
boxplot(mean~grp,data=df_H_CDR3)


summary_CDR3Length_L<-lapply(list_L_CDR3Length,function(x){
  mean<-mean(x$targetColumn)
  median<-median(x$targetColumn)
  q1<-as.numeric(quantile(x$targetColumn,0.1))
  q90<-as.numeric(quantile(x$targetColumn,0.75))
  return(list(mean=mean,median=median,q1=q1,q90=q90))
})
df_L_CDR3<-as.data.frame(unlist(lapply(summary_CDR3Length_L,function(x) x$mean)))
colnames(df_L_CDR3)<-'mean'
df_L_CDR3$grp<-ifelse(row.names(df_L_CDR3) %in% c('Donor_9','Donor_4','Donor_2','Donor_7', 'Donor_5','Donor_6'),'Control','PF')
df_L_CDR3$grp<-as.factor(df_L_CDR3$grp)
boxplot(mean~grp,data=df_L_CDR3)
t.test(mean~grp,data=df_L_CDR3)

#############################################################################################
head(list_H_BestV[[1]])
list_H_BestV_2<-lapply(1:length(list_H_BestV),function(x) {
  list_H_BestV[[x]]$igv<-list_H_BestV[[x]]$targetColumn#
  tab<-as.data.frame(table(list_H_BestV[[x]]$igv))
  colnames(tab)[2]<-names(list_H_BestV)[x];return(tab)})

H_BestV_df<-Reduce(function(x,y) {merge(x,y, by='Var1',all=TRUE)}, list_H_BestV_2);H_BestV_df
row.names(H_BestV_df)<-H_BestV_df$Var1
H_BestV_df<-H_BestV_df[,!colnames(H_BestV_df) %in% 'Var1']
H_BestV_prop<-sweep(H_BestV_df, 2, colSums(H_BestV_df,na.rm=T), FUN = "/");H_BestV_prop

library(tidyr)
H_BestV_prop$IG<-row.names(H_BestV_prop)
H_BestV_prop2<-H_BestV_prop[!(H_BestV_prop$IG) %in% c("IGHV4-28"),]
H_BestV_long <- H_BestV_prop2 %>%
  pivot_longer(
    cols = colnames(H_BestV_prop2)[!colnames(H_BestV_prop2) %in% "IG"], # 
    names_to = "Sample",      
    values_to = "Prop"        
  ) %>% 
  as.data.frame();head(H_BestV_long)
H_BestV_long$Disease<-ifelse(H_BestV_long$Sample %in% c('Donor_9','Donor_4','Donor_2','Donor_7', 'Donor_5','Donor_6'),'Control','PF')
table(H_BestV_long$Disease,H_BestV_long$Sample)

H_BestV_long$Prop[is.na(H_BestV_long$Prop)]<-0
table(H_BestV_long$IG,H_BestV_long$Disease)

(uniIGHV<-unique(H_BestV_long$IG))
P_ighv<-lapply(1:length(uniIGHV),function(x){
  temp<-H_BestV_long[H_BestV_long$IG %in% (uniIGHV[x]),]
  temp$Disease<-as.factor(temp$Disease)
  pt<-t.test(Prop~Disease, data=temp)$p.value
  return(pt)
});names(P_ighv)<-uniIGHV
P_ighv_df<-as.data.frame(unlist(P_ighv))
P_ighv_df$fdr<-p.adjust(P_ighv_df$`unlist(P_ighv)`,'fdr');min(P_ighv_df$fdr)

boxplot(Prop~Disease,data=H_BestV_long[H_BestV_long$IG=='IGHV4-28',])
boxplot(Prop~Disease,data=H_BestV_long[H_BestV_long$IG=='IGHV3-66',])
boxplot(Prop~Disease,data=H_BestV_long[H_BestV_long$IG=='IGHV3-73',])

H_BestV_long[H_BestV_long$IG=='IGHV3-66',]


######################################################################
###########         Public vs Private      ###########################
######################################################################
controlID<-c('Donor_9','Donor_4','Donor_2','Donor_7', 'Donor_5','Donor_6')

H_clon<-Heavy[,grepl("^Fraction.of.UMIs",colnames(Heavy))]
colnames(H_clon)<-sapply(strsplit(colnames(H_clon),split="...",fix=TRUE),function(x) x[2])
H_cln_count<-as.data.frame(apply(H_clon, 1, function(x) sum(!is.na(x))))
colnames(H_cln_count)<-'count'
table(H_cln_count$count)

pdf(sprintf("%s/PublicPrivate_H.pdf",FigOut),width=7.87/2.5,height=7.87/1.2, family="ArialMT")
barplot(table(H_cln_count$count),log='y', col=c('lightblue','tomato','purple'),ylim=c(1,80000))
dev.off()

H_ctrl_clon<-H_clon[,controlID];dim(H_ctrl_clon)
H_pf_clon<-H_clon[,!colnames(H_clon) %in% controlID];dim(H_pf_clon)
Hctrl_cln_count<-as.data.frame(apply(H_ctrl_clon, 1, function(x) sum(!is.na(x))));colnames(Hctrl_cln_count)<-'ctrl';head(Hctrl_cln_count)
Hpf_cln_count<-as.data.frame(apply(H_pf_clon, 1, function(x) sum(!is.na(x))));colnames(Hpf_cln_count)<-'pf';head(Hpf_cln_count)
sum(row.names(Hctrl_cln_count)==row.names(Hpf_cln_count))
H_ctrlpf_count<-cbind(Hctrl_cln_count,Hpf_cln_count);head(H_ctrlpf_count)
H_ctrlpf_count_2<-H_ctrlpf_count[rownames(H_cln_count[H_cln_count$count>1,,drop=F]),]

###############
RAPUP<-c("C-RAPUP", 
         "C-ASJAQ", 
         "C-WVGIW", 
         "C-LGMEZ",
         "C-YQFPO",
         "C-GUCRL",
         "C-RLISB",
         "C-DOESX",
         "C-JWIXG",
         "C-MQPJR",
         "C-SNOGY",
         "C-AAIII",
         "C-DZRZB",
         "C-HNEMT",
         "C-IXWYK",
         "C-WMABI")

H_ctrlpf_count[RAPUP,]
H_clon[RAPUP,]

REDBU<-c("C-REDBU",
         "C-PWCLQ",
         "C-VPWRJ",
         "C-INHNG",
         "C-NVFUL",
         "C-PBHXR",
         "C-PEDTK",
         "C-PXBWO",
         "C-PVJFC",
         "C-RPLEM",
         "C-TUUXV",
         "C-UBAKC",
         "C-VRSUO",
         "C-VTXTK",
         "C-XFJFY",
         "C-XVSQX",
         "C-YLECN",
         "C-EMSJW",
         "C-JUXFI",
         "C-KLAMZ",
         "C-MFRIH",
         "C-OQXDP",
         "C-PNPMK",
         "C-TGAHD",
         "C-VWOGU")

H_ctrlpf_count[REDBU,]
H_clon[REDBU,]
###################



H_ctrlOnly<-sum(H_ctrlpf_count_2$ctrl>0 & H_ctrlpf_count_2$pf==0)
H_both<-sum(H_ctrlpf_count_2$ctrl>0 & H_ctrlpf_count_2$pf>0)
H_pfOnly<-sum(H_ctrlpf_count_2$ctrl==0 & H_ctrlpf_count_2$pf>0)

overay_H<-c(H_ctrlOnly,H_both,H_pfOnly);overay_H
proportions_H <- round(overay_H / sum(overay_H),3)
labels_H <- paste(c('ctrl','both','pf'), " (", proportions_H, "%)", sep="")

pdf(sprintf("%s/PublicPrivate_pie_H.pdf",FigOut),width=7.87,height=7.87, family="ArialMT")
# Create the pie chart
pie(overay_H,
    labels = labels_H,
    main = "Pie Chart with Proportions_H",
    col = rainbow(length(overay_H)))
dev.off()


##################
L_clon<-Light[,grepl("^Fraction.of.UMIs",colnames(Light))];head(L_clon)
colnames(L_clon)<-sapply(strsplit(colnames(L_clon),split="...",fix=TRUE),function(x) x[2]);head(L_clon)
L_cln_count<-as.data.frame(apply(L_clon, 1, function(x) sum(!is.na(x))));head(L_cln_count)
colnames(L_cln_count)<-'count';head(L_cln_count)
table(L_cln_count$count)

pdf(sprintf("%s/PublicPrivate_L.pdf",FigOut),width=7.87/2,height=7.87/1.2, family="ArialMT")
barplot(table(L_cln_count$count),log='y', col=c('lightblue','tomato','purple','green','pink','red'),ylim=c(1,100000))
dev.off()

L_ctrl_clon<-L_clon[,controlID];dim(L_ctrl_clon)
L_pf_clon<-L_clon[,!colnames(L_clon) %in% controlID];dim(L_pf_clon)
Lctrl_cln_count<-as.data.frame(apply(L_ctrl_clon, 1, function(x) sum(!is.na(x))));colnames(Lctrl_cln_count)<-'ctrl'
Lpf_cln_count<-as.data.frame(apply(L_pf_clon, 1, function(x) sum(!is.na(x))));colnames(Lpf_cln_count)<-'pf'
sum(row.names(Lctrl_cln_count)==row.names(Lpf_cln_count))
L_ctrlpf_count<-cbind(Lctrl_cln_count,Lpf_cln_count);head(L_ctrlpf_count)
L_ctrlpf_count_2<-L_ctrlpf_count[rownames(L_cln_count[L_cln_count$count>1,,drop=F]),]

L_ctrlOnly<-sum(L_ctrlpf_count_2$ctrl>0 & L_ctrlpf_count_2$pf==0)
L_both<-sum(L_ctrlpf_count_2$ctrl>0 & L_ctrlpf_count_2$pf>0)
L_pfOnly<-sum(L_ctrlpf_count_2$ctrl==0 & L_ctrlpf_count_2$pf>0)

overay_L<-c(L_ctrlOnly,L_both,L_pfOnly);overay_L
proportions_L <- round(overay_L / sum(overay_L),3)
labels_L <- paste(c('ctrl','both','pf'), " (", proportions_L, "%)", sep="")

pdf(sprintf("%s/PublicPrivate_pie_L.pdf",FigOut),width=7.87,height=7.87, family="ArialMT")
pie(overay_L,
    labels = labels_L,
    main = "Pie Chart with Proportions_L",
    col = rainbow(length(overay_L)))
dev.off()
