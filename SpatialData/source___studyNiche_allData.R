######################
.libPaths(c("/projects/p32215/tools/forR4.4Seu5.1.0",'/home/yyw9094/R/forR4.4Seu5.1.0',"/home/yyw9094/R/x86_64-pc-linux-gnu-library/4.4"))
######################


FindNicheAllData<-function(object, 
                           group.by,
                           assay = "niche", 
                           cluster.name = "niches",
                           buffer=500,
                           neighbors.k = 20, 
                           niches.k = 4,
                           Kiter.max=50){
  allFov<-paste("fov",unique(object$orig.ident),sep="_")
  allCells<-unlist(lapply(1:length(allFov),function(iFov){
    return(Cells(object[[allFov[iFov]]]))
  }))
  if(sum(duplicated(allCells))>0){
    stop('Duplicated Cell ID!!')
  }
  group.labels <- unlist(object[[group.by]][allCells, ])
  groups <- sort(unique(group.labels))
  cell.type.mtx <- matrix(data = 0, nrow = length(allCells), ncol = length(groups))
  rownames(cell.type.mtx) <- allCells
  colnames(cell.type.mtx) <- groups
  cells.idx <- seq_along(allCells)
  group.idx <- match(group.labels, groups)
  cell.type.mtx[cbind(cells.idx, group.idx)] <- 1
  
  maxXcoord=-1e8
  allCoords<-lapply(1:length(allFov),function(iFov){
    tempCoord<-GetTissueCoordinates(object[[allFov[iFov]]], which = "centroids")
    rownames(tempCoord) <- tempCoord[["cell"]]
    tempCoord <- as.matrix(tempCoord[, c("x", "y")])
    
    if(iFov==1){
      FinalCoord=tempCoord
      maxXcoord <<- max(FinalCoord[, "x"])
    }else{
      FinalCoord=tempCoord
      FinalCoord[,'x']<-FinalCoord[,'x']+maxXcoord+buffer
      maxXcoord <<- max(FinalCoord[, "x"])
    };print(maxXcoord)
    return(FinalCoord)
  })
  
  coords<-do.call('rbind',allCoords);head(coords)
  
  if(!sum(row.names(coords)==allCells)){
    stop("CellId not matched!!!")
  }
  
  neighbors <- FindNeighbors(coords, 
                             k.param = neighbors.k, 
                             compute.SNN = FALSE)
  sum.mtx <- as.matrix(neighbors[["nn"]] %*% cell.type.mtx)
  niche.assay <- CreateAssayObject(counts = t(sum.mtx))
  object[[assay]] <- niche.assay
  DefaultAssay(object) <- assay
  object <- ScaleData(object)
  results <- kmeans(x = t(object[[assay]]@scale.data), 
                    centers = niches.k, 
                    iter.max = Kiter.max,
                    nstart = 30)
  object[[cluster.name]] <- results[["cluster"]]
  return(object)
}



