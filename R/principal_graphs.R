# File description --------------------------------------------
#
# This fine contains functions designed to interaqct with VDAOEngine.jar
#
# Only the main function will be exported and made available to the user
#
#



# Configuration Functions --------------------------------------------


#' Define the grammar to produde a circular principal graph
#'
#' This is an internal function and should not be accessed directly by the user.
#'
#' @return A parameter set that can be used to construct a circular principal graph
CircleConfiguration <- function(){

  params <- NULL

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a closed reactangle spanned by first two PCs
  params$initalgorithm = 2

  # universal parameters (should be rarely modified)
  params$epochs[[1]]$id = 1
  params$epochs[[1]]$minimize = 1
  params$epochs[[1]]$eps = 0.01
  params$epochs[[1]]$epsSLAU = 0.01
  params$epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  params$epochs[[1]]$grammartype = 'circle';

  # elasticity parameters
  params$epochs[[1]]$ep = 0.01
  params$epochs[[1]]$rp = 0.1

  # default number of nodes
  params$epochs[[1]]$numiter = 20

  return(params)
}

#' Define the grammar to produde a tree principal graph
#'
#' This is an internal function and should not be accessed directly by the user.
#'
#' @return A parameter set that can be used to construct a tree principal graph
DefaultPrincipalTreeConfiguration <- function(){

  params <- NULL

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a piece of first principal component
  params$initalgorithm = 0

  # universal parameters (should be rarely modified)
  params$epochs[[1]]$id = 1
  params$epochs[[1]]$minimize = 1
  params$epochs[[1]]$eps = 0.01
  params$epochs[[1]]$epsSLAU = 0.01
  params$epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  params$epochs[[1]]$grammartype = 'treeWithPruning';

  # elasticity parameters
  params$epochs[[1]]$ep = 0.01
  params$epochs[[1]]$rp = 0.1

  # default number of nodes
  params$epochs[[1]]$numiter = 20

  return(params)


}

#' Define the grammar to produde a curve (a path) principal graph
#'
#' This is an internal function and should not be accessed directly by the user.
#'
#' @return A parameter set that can be used to construct a curve (a path) principal graph
CurveConfiguration <- function(){

  params <- NULL

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a piece of first principal component
  params$initalgorithm = 0

  # universal parameters (should be rarely modified)
  params$epochs[[1]]$id = 1
  params$epochs[[1]]$minimize = 1
  params$epochs[[1]]$eps = 0.01
  params$epochs[[1]]$epsSLAU = 0.01
  params$epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  params$epochs[[1]]$grammartype = 'curve'

  # elasticity parameters
  params$epochs[[1]]$ep = 0.01
  params$epochs[[1]]$rp = 0.1

  # default number of nodes
  params$epochs[[1]]$numiter = 20

  return(params)
}

#' Define the grammar to produde a robust tree principal graph
#'
#' This is an internal function and should not be accessed directly by the user
#'
#' @return A parameter set that can be used to construct a robust tree principal graph
RobustPrincipalTreeConfiguration <- function(){

  params <- NULL

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a piece of first principal component
  params$initalgorithm = 1

  # universal parameters (should be rarely modified)
  params$epochs[[1]]$id = 1
  params$epochs[[1]]$minimize = 1
  params$epochs[[1]]$eps = 0.01
  params$epochs[[1]]$epsSLAU = 0.01
  params$epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  params$epochs[[1]]$grammartype = 'treeWithPruning';

  # elasticity parameters
  params$epochs[[1]]$ep = 0.01
  params$epochs[[1]]$rp = 0.1

  # default number of nodes
  params$epochs[[1]]$numiter = 30

  # robustness parameters
  params$epochs[[1]]$robust = TRUE
  params$epochs[[1]]$trimradius = 0.1

  return(params)
}

#' Define the grammar to produde a tree principal graph using a two stage approach.
#'
#' This is an internal function and should not be accessed directly by the user.
#' This function is still experimental and not is not working yet.
#'
#' @return A parameter set that can be used to construct a robust tree principal graph
TwoStagesRobustPrincipalTreeConfiguration <- function(){

  params <- list()

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a piece of first principal component
  params$initalgorithm = 0

  ###########
  # Epoch 1 #
  ###########

  epochs <- list()
  epochs[[1]] <- NA

  # universal parameters (should be rarely modified)
  epochs[[1]]$id = 1
  epochs[[1]]$minimize = 1
  epochs[[1]]$eps = 0.01
  epochs[[1]]$epsSLAU = 0.01
  epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  epochs[[1]]$grammartype = 'treeWithPruning';

  # elasticity parameters
  epochs[[1]]$ep = 0.01
  epochs[[1]]$rp = 0.1

  # default number of nodes
  epochs[[1]]$numiter = 30

  # robustness parameters
  epochs[[1]]$robust = FALSE
  epochs[[1]]$trimradius = 0.1

  ###########
  # Epoch 2 #
  ###########

  epochs[[2]] <- NA

  # universal parameters (should be rarely modified)
  epochs[[2]]$id = 1
  epochs[[2]]$minimize = 1
  epochs[[2]]$eps = 0.01
  epochs[[2]]$epsSLAU = 0.01
  epochs[[2]]$numiterSLAU = 1000

  # what grammar to use
  epochs[[2]]$grammartype = 'treeWithPruning';

  # elasticity parameters
  epochs[[2]]$ep = 0.0001
  epochs[[2]]$rp = 0.0001

  # default number of nodes
  epochs[[2]]$numiter = 30

  # robustness parameters
  epochs[[2]]$robust = TRUE
  epochs[[2]]$trimradius = 0.1

  params$epochs <- epochs

  return(params)
}










# Data extraction Functions --------------------------------------------


#' Produce a graph structure from the coputed elastic principal graph details
#'
#' The function depends on a working JVM with the VDAOEngine.jar library in the libpath.
#' This is an internal function and should not be accessed directly by the user
#'
#' @param PrintGraph An elastic principal graph structure as returned from \code{link{computeElasticPrincipalGraph}}
#'
#' @return A graph structure
makeGraph <- function(PrintGraph){

  Graph = .jnew('vdaoengine/analysis/grammars/Graph')

  Utils <- .jnew(class = 'vdaoengine/analysis/grammars/Utils')

  Utils$setNodesForGraph(Graph, .jarray(.jfloat(data.matrix(PrintGraph$Nodes)),dispatch=T))
  if(min(PrintGraph$Edges) == 0){
    Utils$setEdgesForGraph(Graph, .jarray(matrix(as.integer(PrintGraph$Edges), ncol = 2),dispatch=T))
  } else {
    Utils$setEdgesForGraph(Graph, .jarray(matrix(as.integer(PrintGraph$Edges-1), ncol = 2),dispatch=T))
  }

  Graph$defineStarsFromPrimitiveGraphStructure()
  Graph$compileNodesInArrays()

  return(Graph)

}

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
getTaxonMap <- function(Results, Data, UseR = FALSE){

  if(UseR){
    
    CurvePoints <- Results$Nodes
    # rownames(CurvePoints) <- paste("V_", 1:nrow(CurvePoints), sep='')
    
    AllPoints <- rbind(Data, CurvePoints)
    DistPoints <- as.matrix(dist(AllPoints))
    IntDist <- DistPoints[1:nrow(Data),(nrow(Data)+1):(nrow(CurvePoints)+nrow(Data))]
    PointsNodesProjections <- apply(IntDist, 1, which.min)

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
projectPoints <- function(Results, Data, TaxonList=NULL,
                          UseR = TRUE, method = 'Dist'){
  
  if(is.null(TaxonList)){
    print("TaxonList will be computed. Consider doing that separetedly")
    TaxonList <- getTaxonMap(Results = Results, Data = Data, UseR = TRUE)
  }
  
  if(min(Results$Edges) == 0){
    Results$Edges <- Results$Edges + 1
  }
  
  if(UseR){
    
    if(method == 'OnEdges'){
      
      DistsLists <- list()
      
      for(i in 1:nrow(Results$Edges)){
        
        Nd <- Results$Edges[i, ]
        
        C1 <- Results$Nodes[Nd[1],]
        C2 <- Results$Nodes[Nd[2],]
        
        # Select Only points associated with C1 or C2
        
        SelPoints <- c(TaxonList[[Nd[1]]], TaxonList[[Nd[2]]])
        RestDataMap <- Data[SelPoints,]
        
        # plot(RestDataMap[,1:2])
        # points(C1[1], C1[2], pch=20)
        # points(C2[1], C2[2], pch=20)
        # 
        # PTOL <- NULL
        # for(j in seq(from=0, to=1, by=0.01)){
        #   PTOL <- rbind(PTOL, (j)*C1 + (1-j)*C2 )
        # }
        # points(x = PTOL[,1], y = PTOL[,2])
        # 
        # table(apply(as.matrix(dist(rbind(PTOL, RestDataMap)))[1:nrow(PTOL),-c(1:nrow(PTOL))], 2, which.min))
        
        # This is a terrible method. There must be a better way!!!
        
        PCARet <- prcomp(rbind(C1, C2), retx = TRUE)
        
        P1Pos <- PCARet$x["C1","PC1"]
        P2Pos <- PCARet$x["C2","PC1"]
        
        PointPosFull <- t(t(RestDataMap) - PCARet$center) %*% PCARet$rotation
        
        PointPos <- PointPosFull[,"PC1"]
        PrjPoints <- t(t(PointPos %*% t(PCARet$rotation[,1])) + PCARet$center)
        # plot(PrjPoints[,1:2])
        # Nodepoints <- t(t(PCARet$x[,"PC1"] %*% t(PCARet$rotation[,1])) + PCARet$center)
        # points(Nodepoints[,1:2], col='red')
        
        # plot(PrjPoints[,1], PrjPoints[,2], col='red')
        # points(x = PCARet$x[,1], y = PCARet$x[,2])
        
        # arrows(x0 = RestDataMap[,1], y0 = RestDataMap[,2], x1 = PrjPoints[,1], y1 = PrjPoints[,2], length = 0)
        
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
                       SegmentDist = abs(P2Pos-P1Pos),
                       ProjectedCoords = PrjPoints)
        
        DistsLists[[i]] <- ElList
        
      }
      
      PosVector <- matrix(rep(NA, nrow(Data)*(ncol(Data)+1)), nrow = nrow(Data))
      RelPosVector <- matrix(rep(NA, nrow(Data)*3), nrow = nrow(Data))
      SegLen <- NULL
      # plot(Data[,1:2], col="gray", ylim = c(-50, 50))
      # points(Results$Nodes[,1:2], col="black")
      
      for (i in 1:length(DistsLists)) {
        
        SegLen <- c(SegLen, DistsLists[[i]]$SegmentDist)
        
        Inside <- DistsLists[[i]]$PointsProjections > 0 & DistsLists[[i]]$PointsProjections < 1
        
        for (j in 1:length(DistsLists[[i]]$PointsIndices)) {
          if(Inside[j]){
            # The point is inside a different segment. It gets it coordinates
            PosVector[DistsLists[[i]]$PointsIndices[j], ] <- c(1, DistsLists[[i]]$ProjectedCoords[j,])
            RelPosVector[DistsLists[[i]]$PointsIndices[j], ] <-
              c(DistsLists[[i]]$PointsProjections[j], DistsLists[[i]]$Nodes)
            next()
          }
          
          if(!Inside[j] & is.na(PosVector[DistsLists[[i]]$PointsIndices[j], 1])){
            # The point is not part of a different segments.
            # It gets the coordinate of a node
            
            RelPosVector[DistsLists[[i]]$PointsIndices[j], ] <-
              c(DistsLists[[i]]$PointsProjections[j], DistsLists[[i]]$Nodes)
            
            if(DistsLists[[i]]$PointsProjections[j] == 0){
              PosVector[DistsLists[[i]]$PointsIndices[j], ] <- c(0, Results$Nodes[DistsLists[[i]]$Nodes[1],])
            }
            if(DistsLists[[i]]$PointsProjections[j] == 1){
              PosVector[DistsLists[[i]]$PointsIndices[j], ] <- c(0, Results$Nodes[DistsLists[[i]]$Nodes[2],])
            }
            next()
          }
          
        }
        
      }
      
      # for (j in 1:nrow(PosVector)) {
      #   points(x=PosVector[j,2], y=PosVector[j,3], col='green', cex=0.5)
      # }
      
      return(list(PointsOnEdgesCoords = PosVector, EdgeLength = SegLen, PointsOnEdgesDist = RelPosVector))
      
    }
    
    if(method == 'Dist'){
      
      FullMat <- rbind(Results$Nodes, Data)
      DistMat <- as.matrix(dist(FullMat))
      
      RedDistMat <- DistMat[1:nrow(Results$Nodes), -c(1:nrow(Results$Nodes))]
      NodeDistMat <- DistMat[1:nrow(Results$Nodes), c(1:nrow(Results$Nodes))]
      
      if(min(Results$Edges)==0){
        Results$Edges <- Results$Edges + 1
      }
      
      PosVector <- matrix(rep(NA, nrow(Data)*(ncol(Data)+1)), nrow = nrow(Data))
      RelPosVector <- matrix(rep(NA, nrow(Data)*3), nrow = nrow(Data))
      SegLen <- NULL
      
      SortedNodes <- apply(RedDistMat, 2, sort, index.return=TRUE)
      for (i in 1:length(SortedNodes)) {
        BaseNode <- SortedNodes[[i]]$ix[1]
        PossibleNodes <- unique(as.vector(Results$Edges[which(rowSums(Results$Edges == BaseNode)>0),]))
        PossibleNodes <- setdiff(PossibleNodes, BaseNode)
        SecondaryNodeIdx <- min(which(SortedNodes[[i]]$ix %in% PossibleNodes))
        SecondaryNode <- SortedNodes[[i]]$ix[SecondaryNodeIdx]
        
        EdgeId <- which(rowSums(matrix(Results$Edges %in% c(BaseNode, SecondaryNode), ncol = 2)) == 2)
        if(Results$Edges[EdgeId, 1] == BaseNode){
          PointDists <- SortedNodes[[i]]$x[c(1, SecondaryNodeIdx)]
          RelPosVector[i,] <- c(PointDists[1]/sum(PointDists), BaseNode, SecondaryNode)
          PosVector[i, ] <- c(1, RelPosVector[i,1]*Results$Nodes[BaseNode,] +
                                (1-RelPosVector[i,1])*Results$Nodes[SecondaryNode,])
        } else {
          PointDists <- SortedNodes[[i]]$x[c(1, SecondaryNodeIdx)]
          RelPosVector[i,] <- c(PointDists[2]/sum(PointDists), SecondaryNode, BaseNode)
          PosVector[i, ] <- c(1, RelPosVector[i,1]*Results$Nodes[SecondaryNode,] +
                                (1-RelPosVector[i,1])*Results$Nodes[BaseNode,])
        }
        
        SegLen <- c(SegLen, NodeDistMat[BaseNode, SecondaryNode])
        
      }
      
      return(list(PointsOnEdgesCoords = PosVector, EdgeLength = SegLen, PointsOnEdgesDist = RelPosVector))
      
    }
    
    print("Method not supported")
    return(NULL)
       
  } else {
    
    print("Sorry. Not implemented yet ...")
    return(NULL)
    
  }
  
}




# computeElPT --------------------------------------------





#' Compute the the elastic principal graph
#'
#' The function depends on a working JVM with the VDAOEngine.jar library in the libpath.
#' This is an internal function and should not be accessed directly by the user.
#'
#' @param Data The data points to be used to compute the graph
#' @param NumNodes The number of vertices to be used to construct the graph
#' @param Parameters The name of a parameter function
#' @param ... Additional parameter to be procesed by the Java procedure. Currently only EP, RP and TrimRadius are supported
#'
#' @return A list describing the principal elastic graph
computeElPT <- function(Data, NumNodes, Parameters, ...){

  StartTime <- Sys.time()

  EP <- NULL
  RP <- NULL
  TrimRadius <- -1

  AddPar <- list(...)

  if(length(AddPar)>0){

    print("Additional parameters found")
    print(AddPar)

  }

  if(length(AddPar)>0){
    if(any(names(AddPar)=="EP")){
      EP <- AddPar[["EP"]]
      print("EP set from user")
    }

    if(any(names(AddPar)=="RP")){
      RP <- AddPar[["RP"]]
      print("RP set from user")
    }

    if(any(names(AddPar)=="TrimRadius")){
      TrimRadius <- AddPar[["TrimRadius"]]
      print("TrimRadius set from user")
    }

  }


  cat("Configuring engine .")

  Config <- .jnew(class = "vdaoengine/analysis/grammars/ConfigFile")

  Config$algtype <- Parameters$algtype

  tArray <- Config$stretchInitCoeffs
  tArray[1] <- Parameters$initstretch1

  Config$stretchInitCoeffs <- .jfloat(tArray)

  Config$initStrategy <- as.integer(Parameters$initalgorithm)

  cat(".")

  for (i in 1:length(Parameters$epochs)) {

    Epoch <- .jnew(class = "vdaoengine/analysis/elmap/ElmapAlgorithmEpoch")

    Epoch$grammarType <- Parameters$epochs[[i]]$grammartype

    Epoch$EP <- .jfloat(Parameters$epochs[[i]]$ep)

    if(!is.null(EP)){
      if(EP[i]>0){
        Epoch$EP <- .jfloat(EP[i])
      }
    }

    Epoch$RP = .jfloat(Parameters$epochs[[i]]$rp)

    if(!is.null(RP)){
      if(RP[i]>0){
        Epoch$RP <- .jfloat(RP[i])
      }
    }

    Epoch$numberOfIterations <- as.integer(NumNodes[i])
    Epoch$maxNumberOfIterationsForSLAU <- as.integer(Parameters$epochs[[i]]$numiterSLAU)
    Epoch$epsConvergence <- .jfloat(Parameters$epochs[[i]]$eps)
    Epoch$epsConvergenceSLAU <- .jfloat(Parameters$epochs[[i]]$epsSLAU)

    if(!is.null(Parameters$epochs[[1]]$robust)){
      Epoch$robust <- Parameters$epochs[[i]]$robust
    }

    if(!is.null(Parameters$epochs[[1]]$trimradius)){
      Epoch$trimradius <- .jfloat(Parameters$epochs[[i]]$trimradius)
    }

    if(TrimRadius>0){
      Epoch$robust = TRUE
      Epoch$trimradius = TrimRadius[i]
    }

    Config$epochs$add(Epoch)

  }

  cat(".")

  cpg <- .jnew(class = "vdaoengine/analysis/grammars/ComputePrincipalGraph")

  cpg$config <- Config

  cat(".")

  cpg$setDataSetAsMassif(.jarray(.jfloat(data.matrix(Data)),dispatch=T))
  cat(".")

  print("")
  print("Running engine")

  report <- read.delim(text = cpg$compute())

  NodeSize <- cpg$graph$countNumberOfPointsProjected(.jcast(cpg$dataset, new.class = "/java/lang/String" ))

  NodePositions <- cpg$graph$getNodePositions()
  NodePositions <- .jevalArray(NodePositions, simplify = TRUE)

  Edges <- cpg$graph$getEdgeTable()
  Edges <- .jevalArray(Edges, simplify = TRUE)

  EndTime <- Sys.time()


  StructuredReturn <- list()
  StructuredReturn$Nodes <- NodePositions
  StructuredReturn$Edges <- Edges
  StructuredReturn$Report <- report
  StructuredReturn$EndTime <- EndTime
  StructuredReturn$StartTime <- StartTime
  StructuredReturn$NodeSize <- NodeSize


  return(StructuredReturn)
}




# computeElasticPrincipalGraph --------------------------------------------


#' Compute the elastic principal graph
#'
#' @param Data The datapoints to be used. Rows are expected to indicate the points and columes their dimensions
#' @param NumNodes The number of nodes to be used to compute the principal graph
#' @param Method The parametriation associated with a specific graph grammar.
#' @param ReduceMethod An optional string indicating if dimensionality reduction via PCA need to per performed
#' @param ReduceSize Th number of dimensions to retain in case of dimensionality reduction
#' @param ... Additional parameters to be passed to the internal function thqt computes the principal graph
#'
#' @return A list describing the principal elastic graph
#' @export
#'
#' @examples
computeElasticPrincipalGraph <- function(Data, NumNodes, Method = NULL,
                                         ReduceMethod = 'none', ReduceSize = 3, ...) {


  ExtraArgs <- list(...)

  # Reduce dimensions?

  RedDone <- FALSE

  if(ReduceMethod != 'none'){

    print(paste("Performing dimensionality reduction using", ReduceMethod))

    if(ReduceMethod == 'irlba-Lanczos' & ReduceSize == ncol(Data)){

      print(paste("Performing dimensionality reduction using base-svd instead of irlba-Lanczos (RTM)"))

    } else {

      print(paste("Performing dimensionality reduction using", ReduceMethod))

    }

    FunArg <- append(list(Data = Data, Components = ReduceSize, Method = ReduceMethod), ExtraArgs[OkArgs])
    RedData <- do.call(SelectComputePCA, FunArg)

    NewData <- t(t(Data)-RetVal$centers)%*%RetVal$Comp


    if((ReduceMethod != 'base-svd' & ReduceMethod != 'base-cov') & ReduceSize < ncol(Data)){
      print(paste(signif(100*sum(RetVal$ExpVar), 4), "% (estimated) of the original variance vas retained", sep=''))
    } else {
      print(paste(signif(100*sum(RetVal$ExpVar), 4), "% of the original variance vas retained", sep=''))
    }

  } else {

    NewData <- Data
    RedDone <- TRUE

  }


  if(!RedDone){

    print('Invalid dimensionality reduction method selected')
    print('Original data will be used')

  }

  if(is.null(Method)){
    print('Please select the correct parametrization')
  } else {
    Parameters <- do.call(Method, list())
    PrintGraph <- do.call(computeElPT,
                          append(list(Data = NewData, NumNodes = NumNodes, Parameters = Parameters),
                                 ExtraArgs))
    return(PrintGraph)

  }

}

