#' Title
#'
#' @param Results 
#' @param Data 
#' @param Categories 
#' @param Graph 
#' @param Path 
#'
#' @return
#' @export
#'
#' @examples
CheckAndGetPath <- function(Results, Data, Categories, Graph, Path, PathType, Circular) {
  
  # Looking at potential paths
  
  Status <- 'unDone'
  
  if(length(Path) > 1){
    
    for(i in 2:length(Path)){
      if(!igraph::are.connected(Graph, Path[i], Path[i-1])){
        print(paste(Path[i], "--", Path[i-1], "not found in the net"))
        return()
      }
    }
    
    if(Circular){
      if(!igraph::are.connected(Graph, Path[i], Path[i-1])){
        print(paste(Path[length(Path)], "--", Path[1], "not found in the net"))
        return()
      }
    }
    
    if(!Circular){
      if(sum(duplicated(Path))>0){
        print("Duplicated vertices found on path. This is not supported yet")
        return()
      }
    } else {
      if(sum(duplicated(Path[-1]))>0){
        print("Duplicated non-terminal vertices found on path. This is not supported yet")
        return()
      }
    }
    
    PathToUse <- Path
    
  } else {
    
    if(PathType == 'Circular'){
      
      # This assumes that the graph is a closed line
      
      Pattern <- igraph::graph.ring(n = igraph::vcount(Graph), directed = FALSE, mutual = FALSE, circular = TRUE)
      
      PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = Graph, graph2 = Pattern)
      
      if (length(PossiblePaths) == 0) {
        Status <- 'None'
      }
      
      if (length(PossiblePaths) == 1) {
        Status <- 'Single'
      }
      
      if (length(PossiblePaths)> 1) {
        Status <- 'Few'
      }
      
      if (length(PossiblePaths)> 5) {
        Status <- 'Multiple'
      }
      
    }
    
    if(PathType == 'Long.Linear'){
      
      # Looking at all the possible diameters in the graph
      
      DiamLen <- igraph::diameter(graph = Graph, directed = FALSE, unconnected = TRUE)
      
      Pattern <- igraph::graph.ring(n = DiamLen+1, directed = FALSE, mutual = FALSE, circular = FALSE)
      
      PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = Graph, graph2 = Pattern)
      
      if (length(PossiblePaths) == 0) {
        Status <- 'None'
      }
      
      if (length(PossiblePaths) == 1) {
        Status <- 'Single'
      }
      
      if (length(PossiblePaths)> 1) {
        Status <- 'Few'
      }
      
      if (length(PossiblePaths)> 5) {
        Status <- 'Multiple'
      }
      
      if(Circular){
        print("Circular option incompatible with Long.Linear")
      }
      
    }
    
  }
  
  # If multiple path are found ask the user for directions ----------------------------------------------------------
  
  if(Status == 'None'){
    print("No path found")
    return(NULL)
  }
  
  
  if(Status == 'Single'){
    print("A single path was found")
    PathToUse <- PossiblePaths[[1]]$name
  }
  
  
  if(Status == 'Few'){
    print(paste(length(PossiblePaths), "paths found"))
    
    TaxonList <- getTaxonMap(Results = Results, Data = Data, UseR = TRUE)
    
    SelLayOut = 'nicely'
    
    if(PathType == 'Long.Linear'){
      SelLayOut = 'circle_line'
    }
    
    if(PathType == 'Circular'){
      SelLayOut = 'circle'
    }
    
    InfoData <- plotPieNet(Results = Resuls, Graph = Graph,
                           Data = Data, Categories = Categories,
                           TaxonList = TaxonList,
                           NodeSizeMult = 2, ColCat = NULL, LayOut = SelLayOut,
                           DirectionMat = NULL)
    
    legend(x = "center", fill=unique(InfoData$ColInfo[1:3]), legend = unique(Categories))
    
    DONE <- FALSE
    
    while (!DONE) {
      
      print("The following paths have been found:")
      
      for(i in 1:length(PossiblePaths)){
        print(paste("Path", i))
        print(PossiblePaths[[i]])
      }
      
      PathToUseID <- as.integer(readline(prompt="Select the path number that you want to use: "))
      
      if(PathToUseID %in% 1:length(PossiblePaths)){
        print("Accepted")
        DONE <- TRUE
        PathToUse <- PossiblePaths[[PathToUseID]]$name
        print(PathToUse)
      }
      
      if(!DONE){
        print("Invalid Selection.")
      }
    }
  }
  
  
  if(Status == 'Multiple'){
    print(paste(length(PossiblePaths), "paths found"))
    
    TaxonList <- getTaxonMap(Results = Results, Data =  TransfData, UseR = TRUE)
    
    SelLayOut = 'nicely'
    
    if(PathType == 'Long.Linear'){
      SelLayOut = 'circle_line'
    }
    
    if(PathType == 'Circular'){
      SelLayOut = 'circle'
    }
    
    InfoData <- plotPieNet(Results = Results, Graph = Graph,
                           Data = Data, Categories = Categories,
                           TaxonList = TaxonList,
                           NodeSizeMult = 2, ColCat = NULL, LayOut = SelLayOut,
                           DirectionMat = NULL)
    
    legend(x = "center", fill=unique(InfoData$ColInfo[1:3]), legend = unique(CellClass))
    
    DONE <- FALSE
    
    while (!DONE) {
      
      V1 <- readline(prompt="Enter the starting node: ")
      V2 <- readline(prompt="Enter the second node: ")
      
      VertexVect <- unlist(lapply(PossiblePaths, "[[", c(1, 2)), use.names = TRUE)
      VertexMat <- matrix(names(VertexVect), ncol = 2, byrow = TRUE)
      
      PathToUseID <- which(VertexMat[,1] == V1 & VertexMat[,2] == V2)
      
      if(length(PathToUseID) == 1){
        print("Path found:")
        DONE <- TRUE
        PathToUse <- PossiblePaths[[PathToUseID]]$name
        print(PathToUse)
      }
      
      if(!DONE){
        print("Path not found. Please reselect the vertices")
      }
    }
  }
  
  return(PathToUse)
  
}



















#' Title
#'
#' @param ExpressionData 
#' @param CellClass 
#' @param PrinGraph 
#' @param Projections 
#' @param Genes 
#' @param Path 
#'
#' @return
#' @export
#'
#' @examples
GeneExpressiononPath <- function(ExpressionData, TransfData, CellClass = NULL, PrinGraph, Projections, Genes,
                                 InvTransNodes = NULL, Path = 'ask', Net = NULL, PathType = 'Long.Linear',
                                 Circular = FALSE, Plot.Smoother = 'none', Title = "",
                                 Plot = TRUE, Return.Smoother = '', CircExt = .3) {
  
  # Initial checks ----------------------------------------------------------
  
  # Check if cells have a categorization
  
  if(is.null(InvTransNodes) & Plot.Smoother == "Nodes"){
    print("InvTransNodes need to be specifid for Plot.Smoother == Nodes. Using Plot.Smoother = gg")
    Plot.Smoother = "gg"
  }
  
  if(!is.null(InvTransNodes) & Plot.Smoother == 'none'){
    Plot.Smoother = "Nodes"
  }
  
  
  
  if(PathType == 'Long.Linear' & Circular){
    print("Pathtype incompatible with circular options. Using Circular = FALSE")
    Circular <- FALSE
  }
  
  if(CircExt > 1){
    print("Only value of CircExt up to 1 will be considered. Setting CircExt = 1")
    CircExt = 1
  }
  
  if(CircExt <=0){
    print("Only value of CircExt strictly positive will be considered. Setting CircExt = .3")
    CircExt = .3
  }
  
  if(!is.null(CellClass) & length(CellClass) == nrow(ExpressionData)){
    print("Using CellClass to plot classes")
    if(is.null(levels(CellClass))){
      CellClass <- factor(CellClass)
    }
  } else {
    print("CellClass absent or incompatible with data. Classification be ignored")
    CellClass <- factor(rep(1, nrow(ExpressionData)))
  }
  
  FoundGenes <- Genes[Genes %in% colnames(ExpressionData)]
  print(paste(length(FoundGenes), "genes found"))
  
  if(length(FoundGenes) == 0){
    print("Nothing to do")
    return()
  }
  
  if(is.null(Net)){
    print("Constructing graph. Consider doing that separatedly.")
    Net <- ConstructGraph(Results = PrinGraph, DirectionMat = NULL)
  }
  
  
  # Looking for paths if necessary ----------------------------------------------------------
  
  # This section has been moved to a different function
  
  PathToUse <- CheckAndGetPath(Results = PrinGraph, Data = TransfData,
                               Categories = CellClass, Graph = Net, Path = Path,
                               PathType = PathType, Circular = Circular)
  
  
  # Project the points on the path ----------------------------------------------------------
  
  # All preprocessing is done. Now we can look at gene expression over pseudotime
  
  NumericPath <- as.numeric(unlist(lapply(strsplit(PathToUse, "V_"), "[[", 2)))
  
  if(Circular){
    NumericPath <- c(NumericPath, NumericPath[1])
  }
  
  print("Projecting cells on path")
  
  PathProjection <- OrderOnPath(PrinGraph = PrinGraph, Path = NumericPath, PointProjections = Projections)
  
  # PlotOnPath(PathProjection, StageVect)
  
  ProjectedPoints <- which(!is.na(PathProjection$PositionOnPath))
  
  SortedProj <- sort(PathProjection$PositionOnPath[ProjectedPoints], index.return=TRUE)
  
  MatGenesToPlot <- ExpressionData[ProjectedPoints[SortedProj$ix],FoundGenes]

  dim(MatGenesToPlot) <- c(length(ProjectedPoints), length(FoundGenes))
  colnames(MatGenesToPlot) <- FoundGenes
  
  # Prepare data for plotting ----------------------------------------------------------
  
  print("Preparing data")
  
  tictoc::tic()
  
  ggMat <- cbind(rep(SortedProj$x/sum(PathProjection$PathLen), ncol(MatGenesToPlot)),
                 as.vector(MatGenesToPlot),
                 rep(as.character(CellClass[ProjectedPoints[SortedProj$ix]]),
                     ncol(MatGenesToPlot)),
                 rep(colnames(MatGenesToPlot), each = nrow(MatGenesToPlot))
  )
  
  ggMat <- cbind(ggMat, rep("Data", nrow(ggMat)))
  
  if(!is.null(InvTransNodes)){
    NodeMat <- InvTransNodes[NumericPath, FoundGenes]
    
    dim(NodeMat) <- c(length(NumericPath), length(FoundGenes))
    colnames(NodeMat) <- FoundGenes
    
    ggMat.Nodes <- cbind(rep(cumsum(PathProjection$PathLen)/sum(PathProjection$PathLen), ncol(NodeMat)),
                         as.vector(NodeMat),
                         rep(" Nodes", ncol(NodeMat)),
                         rep(colnames(NodeMat), each = nrow(NodeMat))
    )
    
    ggMat.Nodes <- cbind(ggMat.Nodes, rep("Data", nrow(ggMat.Nodes)))
    
    ggMat <- rbind(ggMat, ggMat.Nodes)
    
  }
  
  if(!Circular){
    
    colnames(ggMat) <- c("Pseudo.Time", "Log.Gene.Exp", "Class", "Gene", "Status")
    
    ggMat <- data.frame(ggMat)
    ggMat$Pseudo.Time <- as.numeric(as.character(ggMat$Pseudo.Time))
    ggMat$Log.Gene.Exp <- as.numeric(as.character(ggMat$Log.Gene.Exp))
    
    for(Ref in rev(levels(CellClass))){
      if(Ref %in% levels(ggMat$Class)){
        ggMat$Class <- relevel(ggMat$Class, Ref)
      }
    }
    
    # LM <- lm(ggMat$Log.Gene.Exp ~ ggMat$Pseudo.Time)
    # summary(LM)
    
    if(Plot){
      if(length(unique(ggMat$Class))>1){
        
        p <- ggplot2::ggplot(ggMat[ggMat$Status == "Data" & ggMat$Class != " Nodes",], ggplot2::aes(x = Class, y = Log.Gene.Exp, fill = Class)) +
          ggplot2::geom_boxplot() + ggplot2::facet_wrap(~ Gene) + ggplot2::labs(title = Title, x = "Pseudo time", y = "Gene expression")
        print(p)
        
        p <- ggplot2::ggplot(ggMat, ggplot2::aes(x = Pseudo.Time, y = Log.Gene.Exp)) +
          ggplot2::geom_point(ggplot2::aes(color = Class)) +
          ggplot2::facet_wrap(~ Gene, scales="fixed") +
          ggplot2::labs(x = "Pseudo time", y = "Gene expression") +
          ggplot2::scale_x_continuous(limits = c(0,1))
      } else {
        p <- ggplot2::ggplot(ggMat, ggplot2::aes(x = Pseudo.Time, y = Log.Gene.Exp)) +
          ggplot2::geom_point() +
          ggplot2::facet_wrap(~ Gene, scales="fixed") +
          ggplot2::labs(x = "Pseudo time", y = "Gene expression") +
          ggplot2::scale_x_continuous(limits = c(0,1))
      }
    }
    
  } else {
    
    # We are going to include extra "shadow" points
   
    print("including virtual points")
     
    if(sum(ggMat[, 1]>1-CircExt)>0){
      ggMat_Minus <- ggMat[ggMat[, 1]>1-CircExt, ]
      dim(ggMat_Minus) <- c(sum(ggMat[, 1]>1-CircExt), 5)
      ggMat_Minus[,5] <- "Virtual"
      ggMat_Minus[,1] <- as.numeric(ggMat_Minus[,1]) - 1
      ggMat <- rbind(ggMat, ggMat_Minus)
    }
    
    
    if(sum(ggMat[, 1]<CircExt)>0){
      ggMat_Plus <- ggMat[ggMat[, 1]<CircExt, ]
      dim(ggMat_Plus) <- c(sum(ggMat[, 1]<CircExt), 5)
      ggMat_Plus[,5] <- "Virtual"
      ggMat_Plus[,1] <- as.numeric(ggMat_Plus[,1]) + 1
      ggMat <- rbind(ggMat, ggMat_Plus)
    }

    colnames(ggMat) <- c("Pseudo.Time", "Log.Gene.Exp", "Class", "Gene", "Status")
    
    ggMat <- data.frame(ggMat)
    ggMat$Pseudo.Time <- as.numeric(as.character(ggMat$Pseudo.Time))
    ggMat$Log.Gene.Exp <- as.numeric(as.character(ggMat$Log.Gene.Exp))
    
    for(Ref in rev(levels(CellClass))){
      if(Ref %in% levels(ggMat$Class)){
        ggMat$Class <- relevel(ggMat$Class, Ref)
      }
    }
    
    # LM <- lm(ggMat$Log.Gene.Exp ~ ggMat$Pseudo.Time)
    # summary(LM)
    
    if(Plot){
      if(length(unique(ggMat$Class))>1){
        
        p <- ggplot2::ggplot(ggMat[ggMat$Status == "Data" & ggMat$Class != " Nodes",], ggplot2::aes(x = Class, y = Log.Gene.Exp, fill = Class)) +
          ggplot2::geom_boxplot() + ggplot2::facet_wrap(~ Gene) + ggplot2::labs(title = Title, x = "Pseudo time", y = "Gene expression")
        print(p)
        
        p <- ggplot2::ggplot(ggMat, ggplot2::aes(x = Pseudo.Time, y = Log.Gene.Exp)) +
          ggplot2::geom_point(ggplot2::aes(color = Class, alpha = Status)) +
          ggplot2::facet_wrap(~ Gene, scales="fixed") +
          ggplot2::labs(x = "Pseudo time", y = "Gene expression") +
          ggplot2::scale_alpha_discrete(range = c(1, 0.2)) +
          ggplot2::geom_vline(xintercept = c(0,1), linetype = "dashed") +
          ggplot2::scale_x_continuous(limits = c(-CircExt, 1+CircExt))
      } else {
        p <- ggplot2::ggplot(ggMat, ggplot2::aes(x = Pseudo.Time, y = Log.Gene.Exp)) +
          ggplot2::geom_point(ggplot2::aes(alpha = Status)) +
          ggplot2::facet_wrap(~ Gene, scales="fixed") +
          ggplot2::labs(x = "Pseudo time", y = "Gene expression") +
          ggplot2::scale_alpha_discrete(range = c(1, 0.2)) +
          ggplot2::geom_vline(xintercept = c(0,1), linetype = "dashed") +
          ggplot2::scale_x_continuous(limits = c(-CircExt, 1+CircExt))
      }
    }
    
  }
  
  if(Plot.Smoother == 'gg'){
    p <- p + ggplot2::geom_smooth(data = subset(ggMat, Class != " Nodes"), color = "black", method = "loess")
  }
  
  if(Plot.Smoother == 'Nodes'){
    p <- p + ggplot2::geom_line(data = subset(ggMat, Class == " Nodes"))
  }
  
  tictoc::toc()
  
  if(Plot){
    p <- p + ggplot2::labs(title = Title)
    print("Plotting")
    print(p) 
    
  }
  
  if (ncol(ggMat) == 5) {
    # ggMat <- ggMat[ggMat[,5] == "Real", ]
    ggMat <- ggMat[, -5]
  }
  
  # Return the smoother if necessary ----------------------------------------------------------
  
  if(Return.Smoother == 'lowess'){

    print("Computing lowess smoother")
    
    tictoc::tic()
    
    SpltiData <- split(ggMat[,1:2], ggMat$Gene)
      
    SupportFun <- function(DF) {
      return(lowess(DF$Log.Gene.Exp ~ DF$Pseudo.Time))
    }
    
    Smoother <- lapply(SpltiData, SupportFun)
    
    tictoc::toc()
    
    return(list(Smooth = Smoother, CellsOnPath = PathProjection))
    
  }
  
  if(Return.Smoother == 'loess'){
    
    print("Computing loess smoother")
    
    tictoc::tic()
    
    SpltiData <- split(ggMat[,1:2], ggMat$Gene)
    
    SupportFun <- function(DF) {
      return(loess(DF$Log.Gene.Exp ~ DF$Pseudo.Time))
    }
    
    Smoother <- lapply(SpltiData, SupportFun)
    
    tictoc::toc()
    
    return(list(Smooth = Smoother, CellsOnPath = PathProjection))
    
  }
  
  return(list(Smooth = NULL, CellsOnPath = PathProjection))
  
}







#' Title
#'
#' @param TransfData 
#' @param CellClass 
#' @param PrinGraph 
#' @param Projections 
#' @param Path 
#' @param Net 
#' @param PathType 
#' @param Circular 
#'
#' @return
#' @export
#'
#' @examples
CategoriesOnPath <- function(TransfData, CellClass, PrinGraph, Projections, Path = 'ask',
                             Net = NULL, PathType = 'Long.Linear', Circular = FALSE) {
  
 
  # Initial checks ----------------------------------------------------------
  
  # Check if cells have a categorization

  if(length(unique(CellClass)) <= 1){
    print("Nothing to do")
    return()
  }
  
  if(is.null(Net)){
    print("Constructing graph. Consider doing that separatedly.")
    Net <- ConstructGraph(Results = PrinGraph, DirectionMat = NULL)
  }
  
  
  # Looking for paths if necessary ----------------------------------------------------------
  
  # This section has been moved to a different function
  
  PathToUse <- CheckAndGetPath(Results, Data = TransfData, CellClass, Graph = Net,
                               Path = Path, PathType = PathType, Circular = Circular)
  
  
  # Project the points on the path ----------------------------------------------------------
  
  # All preprocessing is done. Now we can look at gene expression over pseudotime
  
  NumericPath <- as.numeric(unlist(lapply(strsplit(PathToUse, "V_"), "[[", 2)))
  
  if(Circular){
    NumericPath <- c(NumericPath, NumericPath[1])
  }
  
  print("Projecting cells on path")
  
  PathProjection <- OrderOnPath(PrinGraph = PrinGraph, Path = NumericPath, PointProjections = Projections)
  
  # PlotOnPath(PathProjection, StageVect)
  
  ProjectedPoints <- which(!is.na(PathProjection$PositionOnPath))
  
  ggMat <- cbind(PathProjection$PositionOnPath/sum(PathProjection$PathLen), as.character(CellClass))
  colnames(ggMat) <- c("Pseudo.Time", "Class")
  
  ggMat <- data.frame(ggMat)
  ggMat$Pseudo.Time <- as.numeric(as.character(ggMat$Pseudo.Time))
  
  p <- ggplot2::ggplot(ggMat, ggplot2::aes(factor(Class), Pseudo.Time, fill=Class)) +
    ggplot2::geom_boxplot() + ggplot2::coord_flip() +
    ggplot2::labs(y = "Pseudo time", x = "Categories")
  
  print(p)
  
}








#' Title
#'
#' @param ExpressionMatrix 
#' @param NodeOnGenes 
#' @param NodePos 
#' @param CellPos 
#'
#' @return
#' @export
#'
#' @examples
DistanceStatistics <- function(ExpressionMatrix, NodeOnGenes, NodePos, CellPos) {
  
  
  # ExpressionMatrix <- M0Sel$ExpressionData
  # NodeOnGenes <- M0Sel$InvTransNodes
  # NodePos <- cumsum(M0Sel$PathProjection$PathLen)
  # CellPos <- M0Sel$PathProjection$PositionOnPath
  
  NodeOnGenes <- rbind(NodeOnGenes, NodeOnGenes[1,])
  
  GetDistances <- function(i) {
    
    fn <- approxfun(x = NodePos, y = NodeOnGenes[,i], method = 'linear')
    return(fn(CellPos) - ExpressionMatrix[, i])
    
  }
  
  IdxList <- as.list(1:ncol(ExpressionMatrix))
  
  DistList <- lapply(IdxList, GetDistances)
  names(DistList) <- colnames(ExpressionMatrix)
  AbsDistList <- lapply(DistList, abs)
  
  MaxDists <- unlist(lapply(AbsDistList, max))
  MinDists <- unlist(lapply(AbsDistList, min))
  MeanDists <- unlist(lapply(AbsDistList, mean))
  SdDists <- unlist(lapply(AbsDistList, sd))
  IQRDists <- unlist(lapply(AbsDistList, IQR))
  MedianDists <- unlist(lapply(AbsDistList, median))
  
  return(list(Max = MaxDists, Min = MinDists, Median = MedianDists, Sd = SdDists, Mean = MeanDists))
  
}
