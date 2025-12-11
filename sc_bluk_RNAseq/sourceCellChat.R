CellChat_pipeline<-function(sub_CC=seu_H_krt17,type="triMean",trim=NULL){
  ########################################################################################
  ##data input
  #Extract the CellChat input files from a Seurat object
  data.input <- GetAssayData(sub_CC, assay = "RNA", slot = "data") # normalized data matrix
  labels <- sub_CC$CellType
  meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  cellchat <- setIdent(cellchat, ident.use = "group") # set "group" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  (groupSize <- as.numeric(table(cellchat@idents))) # number of cells in each cell group
  ########################################################################################
  ######Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  
  # set the used database in the object
  cellchat@DB <- CellChatDB
  
  ##Preprocessing the expression data for cell-cell communication analysis
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  ########################################################################################
  #Part II: Inference of cell-cell communication network
  ########################################################################################
  #computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.25)
  #Compute the communication probability and infer cellular communication network
  if(type=="triMean"){
    cellchat <- computeCommunProb(cellchat)
  }else if(type=="truncatedMean"){
    cellchat <- computeCommunProb(cellchat,type="truncatedMean",trim=trim)
  }else{
    print("ModifyCodeToAdjustType")
  }
  
  cellchat <- filterCommunication(cellchat, min.cells = 10)# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
 
  #Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellchat);head(df.net);dim(df.net);length(unique(df.net$pathway_name))
  
  # df.net <- subsetCommunication(cellchat,slot.name = "netP");head(df.net);dim(df.net)
  # df.net <- subsetCommunication(cellchat, sources.use ="AberrantEpi", targets.use = 2);head(df.net)
  # df.net <- subsetCommunication(cellchat, signaling = c("FN1", "FGF"));head(df.net)
  # df.net <- subsetCommunication(cellchat, signaling = c("PDGF"));head(df.net)
  # df.net[df.net$source=="AberrantEpi",]
  
  #Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  #Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  uniquePathway<-unique(cellchat@netP$pathways)
  outList<-list(df.net=df.net,
                uniquePathway=uniquePathway,
                cellchat=cellchat)
  return(outList)
}
 
CellChat_pipeline_largeNumber<-function(sub_CC=seu_H_krt17,type="triMean",trim=NULL,min.cells=50){
  ########################################################################################
  ##data input
  #Extract the CellChat input files from a Seurat object
  data.input <- GetAssayData(sub_CC, assay = "RNA", slot = "data") # normalized data matrix
  labels <- sub_CC$CellType
  meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  cellchat <- setIdent(cellchat, ident.use = "group") # set "group" as default cell identity
  levels(cellchat@idents) # show factor levels of the cell labels
  (groupSize <- as.numeric(table(cellchat@idents))) # number of cells in each cell group
  ########################################################################################
  ######Set the ligand-receptor interaction database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  
  # set the used database in the object
  cellchat@DB <- CellChatDB
  
  ##Preprocessing the expression data for cell-cell communication analysis
  cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  ########################################################################################
  #Part II: Inference of cell-cell communication network
  ########################################################################################
  #computeAveExpr(cellchat, features = c("CXCL12","CXCR4"), type =  "truncatedMean", trim = 0.25)
  #Compute the communication probability and infer cellular communication network
  if(type=="triMean"){
    cellchat <- computeCommunProb(cellchat)
  }else if(type=="truncatedMean"){
    cellchat <- computeCommunProb(cellchat,type="truncatedMean",trim=trim)
  }else{
    print("ModifyCodeToAdjustType")
  }
  
  cellchat <- filterCommunication(cellchat, min.cells = min.cells)# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  
  #Extract the inferred cellular communication network as a data frame
  df.net <- subsetCommunication(cellchat);head(df.net);dim(df.net);length(unique(df.net$pathway_name))
  
  # df.net <- subsetCommunication(cellchat,slot.name = "netP");head(df.net);dim(df.net)
  # df.net <- subsetCommunication(cellchat, sources.use ="AberrantEpi", targets.use = 2);head(df.net)
  # df.net <- subsetCommunication(cellchat, signaling = c("FN1", "FGF"));head(df.net)
  # df.net <- subsetCommunication(cellchat, signaling = c("PDGF"));head(df.net)
  # df.net[df.net$source=="AberrantEpi",]
  
  #Infer the cell-cell communication at a signaling pathway level
  cellchat <- computeCommunProbPathway(cellchat)
  #Calculate the aggregated cell-cell communication network
  cellchat <- aggregateNet(cellchat)
  uniquePathway<-unique(cellchat@netP$pathways)
  outList<-list(df.net=df.net,
                uniquePathway=uniquePathway,
                cellchat=cellchat)
  return(outList)
}


cellchat_plot<-function(obj.cellchat=cellchat,outFolder=outFold,objName="H"){
  if (!dir.exists(outFolder)){
    dir.create(outFolder)
  } else {
    print("Dir already exists!")
  }
  groupSize <- as.numeric(table(obj.cellchat@idents))
  
  for(iloop in 1:2){
    if(iloop==1){tiff(sprintf("%s/%s_All.NumOfInteraction.tiff",outFolder,objName),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
    if(iloop==2){pdf(sprintf("%s/%s_All.NumOfInteraction.pdf",outFolder,objName),width=7.1,height=7.1)}
    netVisual_circle(obj.cellchat@net$count, 
                     vertex.weight = groupSize, 
                     weight.scale = T, 
                     label.edge= F, 
                     title.name = "Number of interactions")
    dev.off()
  }
  
  
  for(iloop in 1:2){
    if(iloop==1){tiff(sprintf("%s/%s_All.InteractionStrength.tiff",outFolder,objName),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
    if(iloop==2){pdf(sprintf("%s/%s_All.InteractionStrength.pdf",outFolder,objName),width=7.1,height=7.1)}
    netVisual_circle(obj.cellchat@net$weight, 
                     vertex.weight = groupSize, 
                     weight.scale = T, 
                     label.edge= F, 
                     title.name = "Interaction weights/strength")
    dev.off()
  }
  
  ######Visualize each signaling pathway using Hierarchy plot, Circle plot or Chord diagram
  (uniquePathway<-unique(obj.cellchat@netP$pathways));length(uniquePathway) 
  
  for(ipath in 1:length(uniquePathway)){
    pathways.show <-uniquePathway[ipath]
 
    for(iloop in 1:2){
      ##chord plot
      if(iloop==1){tiff(sprintf("%s/%s_%s.chord.tiff",outFolder,objName,pathways.show),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
      if(iloop==2){pdf(sprintf("%s/%s_%s.chord.pdf",outFolder,objName,pathways.show),width=7.1,height=7.1)}
      netVisual_aggregate(obj.cellchat, signaling = pathways.show, layout = "chord")
      dev.off()
      
      ##circle plot
      if(iloop==1){tiff(sprintf("%s/%s_%s.circle.tiff",outFolder,objName,pathways.show),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
      if(iloop==2){pdf(sprintf("%s/%s_%s.circle.pdf",outFolder,objName,pathways.show),width=7.1,height=7.1)}
      netVisual_aggregate(obj.cellchat, signaling = pathways.show, layout = "circle")
      dev.off()
      
      ##contribution plot
      if(iloop==1){tiff(sprintf("%s/%s_%s.contribution.tiff",outFolder,objName,pathways.show),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
      if(iloop==2){pdf(sprintf("%s/%s_%s.contribution.pdf",outFolder,objName,pathways.show),width=7.1,height=7.1)}
      conp<-netAnalysis_contribution(obj.cellchat, signaling = pathways.show)
      print(conp)
      dev.off()
    }
    
    
    pairLR<- extractEnrichedLR(obj.cellchat, signaling = pathways.show, geneLR.return = FALSE)
    for(iPair in 1:nrow(pairLR)){
      LR.show <- pairLR[iPair,]
      for(iloop in 1:2){
        if(iloop==1){tiff(sprintf("%s/%s_%s.eachPair.%s.tiff",outFolder,objName,pathways.show,LR.show),res=300,width=7.1,height=7.1,compression = "lzw",unit="in")}
        if(iloop==2){pdf(sprintf("%s/%s_%s.eachPair.%s.pdf",outFolder,objName,pathways.show,LR.show),width=7.1,height=7.1)}
        netVisual_individual(obj.cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
        dev.off()
      }
    }
    
  }
}
  

ligand_receptor_source<-function(signaling="PDGF"){
  cc.path.temp<-cc.path_ABC[cc.path_ABC$autocrine=="Para" & cc.path_ABC$pathway_name %in% signaling & (cc.path_ABC$source %in% c("ABCs")),]
  colnames(cc.path.temp)[which(colnames(cc.path.temp)=="GSE")]<-"IDName"
  cc.path.temp$IDName<-sapply(strsplit(gsub("CC_","",cc.path.temp$IDName),split="\\."),function(x) x[1])
  cc.path.join<-dplyr::full_join(cc.path.temp,pCC,by="IDName")
  matrix.wide<-reshape2::dcast(cc.path.join,GSE~interaction_name)
  row.names(matrix.wide)<-matrix.wide$GSE
  matrix.wide<-matrix.wide[,!colnames(matrix.wide) %in% "GSE",drop=FALSE]
  matrix.wide<-matrix.wide[pCC$GSE,]
  matrix.wide$index<-(1:nrow(matrix.wide))
  matrix.wide<-matrix.wide[,!colnames(matrix.wide)=="NA"]
  return(matrix.wide)
}


ligand_receptor_target<-function(signaling="PDGF"){
  cc.path.temp<-cc.path_ABC[cc.path_ABC$autocrine=="Para" & cc.path_ABC$pathway_name %in% signaling & (cc.path_ABC$target %in% c("ABCs")),]
  colnames(cc.path.temp)[which(colnames(cc.path.temp)=="GSE")]<-"IDName"
  cc.path.temp$IDName<-sapply(strsplit(gsub("CC_","",cc.path.temp$IDName),split="\\."),function(x) x[1])
  cc.path.join<-dplyr::full_join(cc.path.temp,pCC,by="IDName")
  matrix.wide<-reshape2::dcast(cc.path.join,GSE~interaction_name)
  row.names(matrix.wide)<-matrix.wide$GSE
  matrix.wide$ID<-row.names(matrix.wide)
  matrix.wide<-matrix.wide[,!colnames(matrix.wide) %in% "GSE",drop=FALSE]
  matrix.wide<-matrix.wide[pCC$GSE,]
  matrix.wide$index<-(1:nrow(matrix.wide))
  matrix.wide<-matrix.wide[,!colnames(matrix.wide)=="NA"]
  return(matrix.wide)
}


calculate_percentage<-function(obj=seuDis,GSE=SelectGSE, genelist=c("PDGFA","PDGFB","PDGFC"),compartment="Epithelial"){
  temp<-subset(obj,Population %in% compartment)
  DefaultAssay(temp)<-"RNA";Idents(temp)<-"CellType"
  
  CountOfExpression<-lapply(GSE,function(iGSE){
    temp_whichGSE<-subset(temp,GSE %in% iGSE)
    newData<-data.frame(t(GetAssayData(temp_whichGSE,slot="counts")[genelist,]))
    newData2<-data.frame(CellType=temp_whichGSE$CellType,sapply(genelist,function(igene){ifelse(newData[,igene],1,0)}))
    
    allCount<-do.call("rbind",lapply(genelist,function(ic){
      temp_ic<-as.data.frame.matrix(table(newData2[,c("CellType",ic)]))
      temp_ic$GSE<-iGSE
      temp_ic$gene<-ic
      temp_ic$CellType<-row.names(temp_ic)
      return(temp_ic)
    }))
    return(allCount)
  })
  data_rbind<-do.call("rbind",CountOfExpression)
  data_rbind$fraction<-data_rbind[,"1"]/(data_rbind[,"1"]+data_rbind[,"0"])
  data_rbind<-data_rbind[order(-data_rbind$fraction),]
  eachCellType<-split(data_rbind,f=data_rbind$gene)
  return(eachCellType)
}
