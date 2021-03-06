# Project points on the graph -----------------------------------------------



#' Project a set of datapoints onto the nodes of principal graph
#'
#' The function depends on a working JVM with the VDAOEngine.jar library in the libpath.
#'
#' @param PrintGraph An elastic principal graph structure as returned from \code{link{computeElasticPrincipalGraph}}
#' @param Data A matrix containing a set of points with the correct number of dimensions. Each row represent a point
#' @param UseR A bolean indicating the use of R (TRUE) of the Java library to perform the projections 
#'
#' @return A list. Each elements of the list represent the vertex of the principal graph and contains a vector reporting
#' the cells (as row number) that are associted with that particular vertex. The vector equal to NA indicates that no
#' cells are associated with that vertex;
#' @export
getTaxonMap <- function(Results, Data, UseR = TRUE){
  
  if(UseR){
    
    CurvePoints <- Results$Nodes
    # rownames(CurvePoints) <- paste("V_", 1:nrow(CurvePoints), sep='')
    
    SelDists <- fields::rdist(Data, CurvePoints)
    
    PointsNodesProjections <- apply(SelDists, 1, which.min)
    
    TaxonMap <- list()
    
    for (i in 1:nrow(Results$Nodes)) {
      Idxs <- which(PointsNodesProjections == i)
      if(length(Idxs)>0){
        TaxonMap[[i]] <- Idxs
      } else {
        TaxonMap[[i]] <- NA
      }
    }
    
    return(TaxonMap)
    
  } else {
    
    Graph <- makeGraph(Results)
    
    NumberOfNodes <- Graph$Nodes$size()
    
    numpoints <- as.integer(nrow(Data))
    coordnum <- as.integer(ncol(Data))
    
    Dataset <- .jnew('vdaoengine/data/VDataSet')
    Dataset$massif <- .jarray(.jfloat(data.matrix(Data)),dispatch=T)
    Dataset$pointCount = numpoints
    Dataset$coordCount = coordnum
    
    Elo = .jnew(class = 'vdaoengine/analysis/grammars/ElasticEnergyOptimization', Dataset, Graph)
    
    Elo$calcTaxons()
    
    TaxonMap <- list()
    
    for (i in 1:NumberOfNodes) {
      Tx <- Elo$taxons$get(as.integer(i-1))
      TaxonSize <- Tx$size()
      
      if(TaxonSize == 0){
        TaxonVect <- NA
      } else {
        TaxonVect <-NULL
        for (j in 1:TaxonSize) {
          TaxonVect[j] <- Tx$get(as.integer(j-1))+1
        }
      }
      
      TaxonMap[[i]] <- TaxonVect
      
    }
    
    return(TaxonMap)
    
  }
  
}





#' Project Points on the edges of the graph
#'
#' @param Results 
#' @param Data 
#' @param TaxonList 
#' @param UseR 
#' @param method 
#'
#' @return
#' @export
#'
#' @examples
projectPoints <- function(Results, Data, TaxonList = NULL, UseR = TRUE, method = 'PCALin', Dims = NULL, Debug = FALSE){
  
  if(is.null(Dims)){
    Dims <- ncol(Data)
  }
  
  if(is.null(TaxonList)){
    print("TaxonList will be computed. Consider doing that separetedly")
    TaxonList <- getTaxonMap(Results = Results, Data = Data, UseR = TRUE)
  }
  
  if(min(Results$Edges) == 0){
    Results$Edges <- Results$Edges + 1
  }
  
  if(UseR){
    
    if(method == 'PCALin'){
      
      DistsLists <- list()
      
      for(i in 1:nrow(Results$Edges)){
        
        # For each edge
        
        # Get the nodes connecting the edges and reduce their dimension using Dim
        
        Nd <- Results$Edges[i, ]
        C1 <- Results$Nodes[Nd[1],1:Dims]
        C2 <- Results$Nodes[Nd[2],1:Dims]
        
        # Select Only points associated with either C1 or C2
        # To do that I look at the association between points and nodes
        
        SelPoints <- c(TaxonList[[Nd[1]]], TaxonList[[Nd[2]]])
        SelPoints <- SelPoints[!is.na(SelPoints)]
        
        # This is probaly a terrible method, but it's rather quick and straightforward
        # Compute the PCA of the two nodes (All the variance will be explained by the 1st component)
        
        PCARet <- prcomp(rbind(C1, C2), retx = TRUE, center = TRUE, scale. = FALSE)
        
        P1Pos <- PCARet$x["C1","PC1"]
        P2Pos <- PCARet$x["C2","PC1"]
        
        if(P1Pos < P2Pos){
          SLen <- P2Pos-P1Pos
        } else {
          SLen <- P1Pos-P2Pos
        }
        
        if(length(SelPoints)>0){
          RestDataMap <- Data[SelPoints,1:Dims]
          
          # Use the rotation matrix obtained by the PCA to project the nodes.
          
          if(length(SelPoints)>1){
            PointPosFull <- t(t(RestDataMap) - PCARet$center) %*% PCARet$rotation
          } else {
            PointPosFull <- t(RestDataMap - PCARet$center) %*% PCARet$rotation
          }
          
          # The 1st component will allow me to hortogonally project the points on the lines
          
          # PointPos <- PointPosFull[,"PC1"]
          
          PointPos <- PointPosFull
          PointPos[, -1] <- 0
          
          # Get the position on the line by reversing the PCA
          
          # PrjPoints <- t(t(PointPos %*% t(PCARet$rotation[,1])))
          
          # t(t(pca$x %*% t(pca$rotation)) + pca$center)
          
          PrjPoints <- t(t(PointPos %*% t(PCARet$rotation)) + PCARet$center)
          
          # Plotting stuff to debug
          
          # plot(PrjPoints[,1:2])
          # Nodepoints <- t(t(PCARet$x[,"PC1"] %*% t(PCARet$rotation[,1])) + PCARet$center)
          # points(Nodepoints[,1:2], col='red')
          
          # scatterplot3d::scatterplot3d(rbind(PrjPoints[,1:3], C1[1:3], C2[1:3]), color = c(rep("black", nrow(PrjPoints)), "red", "red"))
          # 
          # rgl::plot3d(rbind(PrjPoints[,1:3], C1[1:3], C2[1:3]), col = c(rep("black", nrow(PrjPoints)), "red", "red"))
          
          # plot(PrjPoints[,1], PrjPoints[,2], col='red')
          # points(x = PCARet$x[,1], y = PCARet$x[,2])
          
          # arrows(x0 = RestDataMap[,1], y0 = RestDataMap[,2], x1 = PrjPoints[,1], y1 = PrjPoints[,2], length = 0)
          
          # Due to the high dimensionality of the system some poits will project outside of the
          # segment joining the two nodes. In this case the points are associated with the nearest node.
          
          if(Debug){
            print(paste(P1Pos, P2Pos))
            print(PointPos)
          }
          
          if(P1Pos < P2Pos){
            PointPos <- (PointPos - P1Pos) / (P2Pos-P1Pos)
            PointPos[PointPos < 0] <- 0
            PointPos[PointPos > 1] <- 1
          } else {
            PointPos <- (PointPos - P2Pos) / (P1Pos-P2Pos)
            PointPos[PointPos < 0] <- 0
            PointPos[PointPos > 1] <- 1
            PointPos <- 1 - PointPos
          }
          
          ElList <- list(Nodes = Nd,
                         PointsProjections = PointPos,
                         PointsIndices = SelPoints,
                         SegmentDist = SLen,
                         ProjectedCoords = PrjPoints)
          
        } else {
          
          ElList <- list(Nodes = Nd,
                         PointsProjections = NULL,
                         PointsIndices = NULL,
                         SegmentDist = SLen,
                         ProjectedCoords = NULL)
          
        }
        
        DistsLists[[i]] <- ElList
        
      }
      
      # Lets put all together
      
      # PosVect will contain the projections of the points on the edges
      PosVector <- matrix(rep(NA, nrow(Data)*(Dims)), nrow = nrow(Data))
      
      # RelPosVector will contain the relative position on the edges, the nodes that control the edge,
      # the length of the distance between the projection and the point, and the point segment
      RelPosVector <- matrix(rep(NA, nrow(Data)*5), nrow = nrow(Data))
      
      # SegLen will contain the length of the segment. This is restricted to the number of dimensions selected by Dims.
      # Therefore the length will change depending on Dim
      SegLen <- NULL
      
      # onEdge keeps track if the point has been projected on a node (0) or edge (1)
      onEdge <- rep(NA, nrow(Data))
      
      # Plotting for debugging purposes
      # plot(Data[,1:2], col="gray", ylim = c(-50, 50))
      # points(Results$Nodes[,1:2], col="black")
      
      for (i in 1:length(DistsLists)) {
        
        SegLen <- c(SegLen, DistsLists[[i]]$SegmentDist)
        
        if(is.null(DistsLists[[i]]$PointsProjections)){
          next()
        }
        
        Inside <- DistsLists[[i]]$PointsProjections > 0 & DistsLists[[i]]$PointsProjections < 1
        
        if(Debug){
          print(Inside)
        }
        
        for (j in 1:length(DistsLists[[i]]$PointsIndices)) {
          if(Inside[j]){
            # The point is inside a segment. I will compare distances anyway
            
            if(Debug){
              print(paste(j, Inside[j]))
            }
            
            # Has the point been assigned before ?
            
            Assign <- TRUE
            
            NewDist <- sqrt(
              sum(
                (DistsLists[[i]]$ProjectedCoords[j,] - Data[DistsLists[[i]]$PointsIndices[j],1:Dims])^2
              )
            )
            
            if(!is.na(onEdge[DistsLists[[i]]$PointsIndices[j]])){
              # There is a previous assignement. I need to check distances
              # First of all I'm going to obtain the old distance
              OldDist <- sqrt(
                sum(
                  (PosVector[DistsLists[[i]]$PointsIndices[j], ] - Data[DistsLists[[i]]$PointsIndices[j],1:Dims])^2
                )
              )
              
              if(NewDist>OldDist){
                # New distance is larger I will not reassign the point
                Assign <- FALSE
              } else{
                if(Debug){
                  print(paste("Updating projection of point", DistsLists[[i]]$PointsIndices[j]))
                }
                NewDist <- OldDist
              }
              
            }
            
            if(Assign){
              if(Debug){
                print(paste("Point", DistsLists[[i]]$PointsIndices[j], "assigned to edge", paste(DistsLists[[i]]$Nodes, collapse=' ')))
              }
              onEdge[DistsLists[[i]]$PointsIndices[j]] <- 1
              PosVector[DistsLists[[i]]$PointsIndices[j], ] <- DistsLists[[i]]$ProjectedCoords[j,]
              RelPosVector[DistsLists[[i]]$PointsIndices[j], ] <-
                c(DistsLists[[i]]$PointsProjections[j],
                  DistsLists[[i]]$Nodes,
                  NewDist,
                  SegLen[length(SegLen)])
            }
            
            # I'm done with this point. On to the next one
            next()
          }
          
          if(!Inside[j]){
            # The point has not been projected on a segment.
            # It gets the coordinate of a node
            
            if(Debug){
              print(paste("Point", DistsLists[[i]]$PointsIndices[j], "assigned to node", paste(DistsLists[[i]]$Nodes, collapse=' ')))
            }
            
            # Has the point been assigned before ?
            
            Assign <- TRUE
            
            NewDist <- sqrt(
              sum(
                (DistsLists[[i]]$ProjectedCoords[j,] - Data[DistsLists[[i]]$PointsIndices[j],1:Dims])^2
              )
            )
            
            if(!is.na(onEdge[DistsLists[[i]]$PointsIndices[j]])){
              # There is a previous assignement. I need to check distances
              # First of all I'm going to obtain the old distance
              OldDist <- sqrt(
                sum(
                  (PosVector[DistsLists[[i]]$PointsIndices[j], ] - Data[DistsLists[[i]]$PointsIndices[j],1:Dims])^2
                )
              )
              
              if(NewDist>OldDist){
                # New distance is larger I will not reassign the point
                Assign <- FALSE
              } else{
                if(Debug){
                  print(paste("Updating projection of point", DistsLists[[i]]$PointsIndices[j]))
                }
                NewDist <- OldDist
              }
              
            }
            
            if(DistsLists[[i]]$PointsProjections[j] == 0){
              if(Debug){
                print("Node 0")
              }
              onEdge[DistsLists[[i]]$PointsIndices[j]] <- 0
              PosVector[DistsLists[[i]]$PointsIndices[j], ] <- Results$Nodes[DistsLists[[i]]$Nodes[1],1:Dims]
            }
            if(DistsLists[[i]]$PointsProjections[j] == 1){
              if(Debug){
                print("Node 1")
              }
              onEdge[DistsLists[[i]]$PointsIndices[j]] <- 0
              PosVector[DistsLists[[i]]$PointsIndices[j], ] <- Results$Nodes[DistsLists[[i]]$Nodes[2],1:Dims]
            }
            
            
            RelPosVector[DistsLists[[i]]$PointsIndices[j], ] <-
              c(DistsLists[[i]]$PointsProjections[j],
                DistsLists[[i]]$Nodes,
                NewDist,
                SegLen[length(SegLen)])
            
            # I'm done. On to the next point (not really needed, just for simmetry)
            next()
          }
          
        }
      }
      
      # for (j in 1:nrow(PosVector)) {
      #   points(x=PosVector[j,2], y=PosVector[j,3], col='green', cex=0.5)
      # }
      
      return(list(Dims = Dims, PointsOnEdgesCoords = PosVector, OnEdge = onEdge, EdgeLength = SegLen, PointsOnEdgesDist = RelPosVector))
      
    }
    
    if(method == 'Dist'){
      
      print("This method has not been implemented yet!")
      
      # FullMat <- rbind(Results$Nodes[,1:Dims], Data[,1:Dims])
      # DistMat <- as.matrix(dist(FullMat))
      # 
      # RedDistMat <- DistMat[1:nrow(Results$Nodes), -c(1:nrow(Results$Nodes))]
      # NodeDistMat <- DistMat[1:nrow(Results$Nodes), c(1:nrow(Results$Nodes))]
      # 
      # if(min(Results$Edges)==0){
      #   Results$Edges <- Results$Edges + 1
      # }
      # 
      # PosVector <- matrix(rep(NA, nrow(Data)*(Dims+1)), nrow = nrow(Data))
      # RelPosVector <- matrix(rep(NA, nrow(Data)*3), nrow = nrow(Data))
      # SegLen <- NULL
      # 
      # SortedNodes <- apply(RedDistMat, 2, sort, index.return=TRUE)
      # for (i in 1:length(SortedNodes)) {
      #   BaseNode <- SortedNodes[[i]]$ix[1]
      #   PossibleNodes <- unique(as.vector(Results$Edges[which(rowSums(Results$Edges == BaseNode)>0),]))
      #   PossibleNodes <- setdiff(PossibleNodes, BaseNode)
      #   SecondaryNodeIdx <- min(which(SortedNodes[[i]]$ix %in% PossibleNodes))
      #   SecondaryNode <- SortedNodes[[i]]$ix[SecondaryNodeIdx]
      #   
      #   EdgeId <- which(rowSums(matrix(Results$Edges %in% c(BaseNode, SecondaryNode), ncol = 2)) == 2)
      #   if(Results$Edges[EdgeId, 1] == BaseNode){
      #     PointDists <- SortedNodes[[i]]$x[c(1, SecondaryNodeIdx)]
      #     RelPosVector[i,] <- c(PointDists[1]/sum(PointDists), BaseNode, SecondaryNode)
      #     PosVector[i, ] <- c(1, RelPosVector[i,1]*Results$Nodes[BaseNode,1:Dims] +
      #                           (1-RelPosVector[i,1])*Results$Nodes[SecondaryNode,1:Dims])
      #   } else {
      #     PointDists <- SortedNodes[[i]]$x[c(1, SecondaryNodeIdx)]
      #     RelPosVector[i,] <- c(PointDists[2]/sum(PointDists), SecondaryNode, BaseNode)
      #     PosVector[i, ] <- c(1, RelPosVector[i,1]*Results$Nodes[SecondaryNode,1:Dims] +
      #                           (1-RelPosVector[i,1])*Results$Nodes[BaseNode,1:Dims])
      #   }
      #   
      #   SegLen <- c(SegLen, NodeDistMat[BaseNode, SecondaryNode])
      #   
      # }
      # 
      # return(list(PointsOnEdgesCoords = PosVector, EdgeLength = SegLen, PointsOnEdgesDist = RelPosVector))
      
    }
    
    print("Method not supported")
    return(NULL)
    
  } else {
    
    print("Sorry. Not implemented yet ...")
    return(NULL)
    
  }
  
}







# Order points on a path -----------------------------------------------

#' Title
#'
#' @param PrinGraph 
#' @param Path 
#' @param PointProjections 
#'
#' @return
#' @export
#'
#' @examples
OrderOnPath <- function(PrinGraph, Path, PointProjections){
  
  # Check that the path start from 1
  
  if(min(PrinGraph$Edges)==0){
    PrinGraph$Edges <- PrinGraph$Edges + 1
  }
  
  # Trasverse the path
  
  Basedist <- 0
  PathLen <- 0
  
  PointPos <- rep(NA, length(PointProjections$OnEdge))
  PointDist <- rep(NA, length(PointProjections$OnEdge))
  Indices <- rep(NA, length(PointProjections$OnEdge))
  
  for(i in 2:length(Path)){
    
    EdgId <- which(unlist(lapply(apply(PrinGraph$Edges, 1, intersect, Path[c(i-1, i)]), length)) == 2)
    
    if(length(EdgId) != 1){
      stop("Path not found in the graph")
    }
    
    IdPtOnEdge <- which(unlist(lapply(apply(PointProjections$PointsOnEdgesDist[,2:3], 1, intersect, Path[c(i-1, i)]), length))==2)
    
    print(paste(length(IdPtOnEdge), "points found on edge"))
    
    if(length(IdPtOnEdge)>1){
      SubInfo <- PointProjections$PointsOnEdgesDist[IdPtOnEdge,]
      
      if(all(Path[c(i-1, i)] == PrinGraph$Edges[EdgId,])){
        # Path adn edge in the principla graph have the same direction 
        ToOrder <- SubInfo[,1]*SubInfo[,5]
      } else {
        # Path adn edge in the principla graph have the opposite direction 
        ToOrder <- (1-SubInfo[,1])*SubInfo[,5]
      }
      
      SortData <- sort(ToOrder, index.return=TRUE)
      
      PointPos[IdPtOnEdge[SortData$ix]] <- SortData$x + sum(Basedist)
      PointDist[IdPtOnEdge[SortData$ix]] <- SubInfo[SortData$ix,4]
    }
    
    if(length(IdPtOnEdge)==1){
      SubInfo <- PointProjections$PointsOnEdgesDist[IdPtOnEdge,]
      
      if(all(Path[c(i-1, i)] == PrinGraph$Edges[EdgId,])){
        # Path adn edge in the principla graph have the same direction 
        PointPos[IdPtOnEdge] <- SubInfo[1]*SubInfo[5] + sum(Basedist)
      } else {
        # Path adn edge in the principla graph have the opposite direction 
        PointPos[IdPtOnEdge] <- (1-SubInfo[1])*SubInfo[5] + sum(Basedist)
      }
      
      PointDist[IdPtOnEdge] <- SubInfo[4]
    }

    Basedist <- c(Basedist, PointProjections$EdgeLength[EdgId])
    
  }
  
  return(list(PositionOnPath = PointPos, DistanceFromPath = PointDist, PathLen = Basedist))

}

