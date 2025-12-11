######################
.libPaths(c("/projects/p32215/tools/forR4.4Seu5.1.0",'/home/yyw9094/R/forR4.4Seu5.1.0',"/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4"))
######################

path<-"/projects/b1042/Yuanqing/SG/DataDown/GSE150910";list.files(path)
dat1<-read.csv(sprintf("%s/GSE150910_gene-level_count_file.csv",path),row.names=1)
##############################
###source("./bulkRNAseq_sourceCode.R")
pdata<-data.frame(ID=colnames(dat1),
                  type=sapply(strsplit(colnames(dat1),split="\\_"),function(x) x[1]))
pdata$type <- factor(pdata$type, levels=c("control","chp","ipf"))
dat_input<-dat1[,paste(pdata$ID)]
dat_input<-dat_input[rowSums(dat_input)>0,]

library(edgeR)
edgeInput1<-DGEList(dat_input,group=as.factor(pdata$type))
edgeInput2 <- edgeInput1[rowSums(1e+06*edgeInput1$counts/expandAsMatrix(edgeInput1$samples$lib.size,dim(edgeInput1))>=1)>=ncol(dat_input)*0.1, ]
edgeInput2 <- calcNormFactors(edgeInput2)
#############################################################################################################
#############################################################################################################
design.mat <- model.matrix(~pdata$type)
colnames(design.mat) <- sapply(strsplit(colnames(design.mat),split="type"),function(x) x[2])
colnames(design.mat)[c(1)]<- c("(Intercept)")
rownames(design.mat) <- colnames(edgeInput2)
#############################################################################################################
#############################################################################################################
edgeInput3 <- estimateDisp(edgeInput2, design.mat, robust=TRUE) ##
edgeInput3$common.dispersion;#plotBCV(edgeInput3)
fit <- glmQLFit(edgeInput3, design.mat, robust=TRUE);#plotQLDisp(fit)
#############################################################################################################
chpvsCtrl <- glmQLFTest(fit, coef=2)
chpvsCtrl$table$FDR<-p.adjust(chpvsCtrl$table$PValue,"BH")
ipfvsCtrl <- glmQLFTest(fit, coef=3)
ipfvsCtrl$table$FDR<-p.adjust(ipfvsCtrl$table$PValue,"BH")
ipfvschp <- glmQLFTest(fit, contrast=c(0,-1,1))
ipfvschp$table$FDR<-p.adjust(ipfvschp$table$PValue,"BH")
###########
cpm_data<-cpm(edgeInput2,prior.count=1,log=TRUE)

AllOut<-list(dat1=dat1,
             fit=fit,
             chpvsCtrl=chpvsCtrl,
             ipfvsCtrl=ipfvsCtrl,
             ipfvschp=ipfvschp,
             cpm=cpm_data)

saveRDS(AllOut, file="GSE150910_DE_Result.RDS")


##############################################################################################################
##############################################################################################################
##############################################################################################################
AllOut<-readRDS(file="GSE150910_DE_Result.RDS")
dat1=AllOut$dat1
fit=AllOut$fit
chpvsCtrl=AllOut$chpvsCtrl
ipfvsCtrl=AllOut$ipfvsCtrl
ipfvschp=AllOut$ipfvschp
cpm=AllOut$cpm_data

pltTab<-ipfvsCtrl$table;sum(pltTab$FDR<0.05)
pltTab$col<-ifelse(pltTab$FDR<0.05 & abs(pltTab$logFC)>log2(2),'red','black');sum(pltTab$col=='red')
genIn<-c("CXCL13","FDCSP","LY6D","IGLL5","MZB1","FCRLA","FCRL5",'CD79A','MS4A1','TNFRSF13C')
geneLab=pltTab[genIn,]

pdf("/projects/b1042/Yuanqing/MA_bulkRNAseq.pdf", width=8,height=7)
plot(ipfvsCtrl$table$logCPM, ipfvsCtrl$table$logFC,xlab="Log CPM",ylab="Log FC",col=pltTab$col)
abline(h=0,lty=2,col='grey')
text(geneLab$logCPM,geneLab$logFC,labels=row.names(geneLab),adj=0)
dev.off()


pltTabHP<-chpvsCtrl$table;sum(pltTabHP$FDR<0.05)
pltTabHP<-pltTabHP[order(-pltTabHP$logFC),]
pltTabHP$col<-ifelse(pltTabHP$FDR<0.05 & abs(pltTabHP$logFC)>log2(2),'red','black')
genInHP<-c("CXCL13","FDCSP","LY6D","IGLL5","MZB1","FCRLA","FCRL5",'CD79A','MS4A1','TNFRSF13C')
geneLabHP=pltTabHP[genInHP,]

pdf("/projects/b1042/Yuanqing/MA_bulkRNAseq_cHP.pdf", width=8,height=7)
plot(chpvsCtrl$table$logCPM, chpvsCtrl$table$logFC,xlab="Log CPM",ylab="Log FC",col=pltTabHP$col)
abline(h=0,lty=2,col='grey')
text(geneLabHP$logCPM,geneLabHP$logFC,labels=row.names(geneLabHP),adj=0)
dev.off()

pltTabIPFHP<-ipfvschp$table;sum(pltTabIPFHP$FDR<0.05)
pltTabIPFHP$col<-ifelse(pltTabIPFHP$FDR<0.05 & abs(pltTabIPFHP$logFC)>log2(2),'red','black');sum(pltTabIPFHP$col=='red')

ipfUp<-pltTab[pltTab$col=='red' & pltTab$logFC>0,];dim(ipfUp)
ipfDn<-pltTab[pltTab$col=='red' & pltTab$logFC<0,];dim(ipfDn)
chpUp<-pltTabHP[pltTabHP$col=='red' & pltTabHP$logFC>0,];dim(chpUp)
chpDn<-pltTabHP[pltTabHP$col=='red' & pltTabHP$logFC<0,];dim(chpDn)


library(nVennR)
gene_list<-list(ipfUp=row.names(ipfUp),
                ipfDn=row.names(ipfDn),
                chpUp=row.names(chpUp),
                chpDn=row.names(chpDn))

myNV <- plotVenn(gene_list, 
                 outFile = "/projects/b1042/Yuanqing/bulkIPFcHP.svg", 
                 systemShow = TRUE)

bothUp<-intersect(row.names(ipfUp),row.names(chpUp))

library(clusterProfiler)
library(org.Hs.eg.db)

gene_df <- bitr(
  bothUp,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db)

ego_bp <- enrichGO(
  gene          = gene_df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 0.3,
  readable      = TRUE  # 
);dim(ego_bp)

bpro<-ego_bp[ego_bp$Description %in% c('B cell proliferation'),]



##############################################################################################################
##############################################################################################################
##############################################################################################################

AllOut<-readRDS(file="GSE150910_DE_Result.RDS")
dat1=AllOut$dat1
fit=AllOut$fit
chpvsCtrl=AllOut$chpvsCtrl
ipfvsCtrl=AllOut$ipfvsCtrl
ipfvschp=AllOut$ipfvschp
cpm=AllOut$cpm_data

ipfTab<-ipfvsCtrl$table;head(ipfTab)
ipfTab$pos<-ifelse(ipfTab$logFC>1 & ipfTab$FDR<0.05,'Pos','UN')
ipfTab$pos<-ifelse(ipfTab$logFC<(-1) & ipfTab$FDR<0.05,'Neg',ipfTab$pos);table(ipfTab$pos)
ipfTab$stat_gsea<-zscoreT(sign(ipfTab$logFC)*sqrt(ipfTab$F),df=ipfvsCtrl$df.total)

library(msigdbr)
library(fgsea)
library(ggplot2)
library(data.table)
require(gridExtra)

ranks <- ipfTab$stat_gsea;head(ranks)
names(ranks) <- row.names(ipfTab)
ranks<-ranks[order(-ranks)];head(ranks)

# Running fgsea
stats <- ranks
geneSets <- msigdbr(species = "Homo sapiens", category = 'C5', subcategory = "GO:BP")
geneSets <- geneSets[geneSets$gene_symbol %in% names(ranks),]
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)

set.seed(135)
eaRes <- fgsea(pathways = m_list, stats = stats, nperm = 5e4, minSize = 15)
eaRes <-eaRes[order(eaRes$padj),]
fwrite(eaRes,file="/projects/b1042/Yuanqing/RNAseq_GOBP_GSEA.csv")
fwrite(eaRes[eaRes$padj<0.05,],file="/projects/b1042/Yuanqing/RNAseq_GOBP_GSEA_Sig.csv")

selBPPath<-c("GOBP_B_CELL_ACTIVATION",
             'GOBP_B_CELL_PROLIFERATION',
             'GOBP_HUMORAL_IMMUNE_RESPONSE',
             'GOBP_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN')

fourBP<-lapply(1:length(selBPPath),function(x){
  return(plotEnrichment(m_list[[selBPPath[x]]], stats) +
          labs(title=selBPPath[x]))
})

pdf("/projects/b1042/Yuanqing/bulk_4BP.pdf",width=7.87,height=7.87/1.2, family="ArialMT")
do.call(grid.arrange, c(fourBP, ncol=2,nrow=2))
dev.off()
