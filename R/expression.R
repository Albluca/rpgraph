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
                                 Path = 'ask', Net = NULL, PathType = 'Long.Linear', Circular = FALSE) {
  
  # Check if cells have a categorization
  
  if(!is.null(CellClass) & length(CellClass) == nrow(ExpressionData)){
    print("Using CellClass to plot classes")
    if(is.null(levels(CellClass))){
      CellClass <- factor(CellClass)
    }
  } else {
    print("CellClass absent or incompatible with data. Classification be ignored")
    CellClass <- factor(rep(1, nrow(ExpressionData)))
  }
  
  FoundGenes <- Genes %in% colnames(ExpressionData)
  print(paste(sum(FoundGenes), "genes found"))
  
  if(sum(FoundGenes) == 0){
    print("Nothing to do")
    return()
  }
  
  # Looking at potential paths
  
  if(is.null(Net)){
    print("Constructing graph. Consider doing that separatedly.")
    Net <- ConstructGraph(Results = PrinGraph, DirectionMat = NULL)
  }
  
  
  
  Status <- 'unDone'
  
  if(length(Path) > 1){
    
    for(i in 2:length(Path)){
      if(!igraph::are.connected(Net, Path[i], Path[i-1])){
        print(paste(Path[i], "--", Path[i-1], "not found in the net"))
        return()
      }
    }
    
    if(Circular){
      if(!igraph::are.connected(Net, Path[i], Path[i-1])){
        print(paste(Path[length(Path)], "--", Path[1], "not found in the net"))
        return()
      }
    }
    
    if(!Circular){
      if(sum(duplicated(Path))>0){
        print("Duplicqted vertex found on path. This is not supported yet")
        return()
      }
    } else {
      if(sum(duplicated(Path[-1]))>0){
        print("Duplicated non-terminal vertex found on path. This is not supported yet")
        return()
      }
    }
    
  } else {
    
    if(PathType == 'Circular'){
      
      # This assules that the graph is a closed line
      
      Pattern <- igraph::graph.ring(n = igraph::vcount(Net), directed = FALSE, mutual = FALSE, circular = TRUE)
      
      PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = Net, graph2 = Pattern)
      
      if (length(PossiblePaths) == 0) {
        Status <- 'None'
      }
      
      if (length(PossiblePaths) == 1) {
        Status <- 'Single'
      }
      
      if (length(PossiblePaths)> 1) {
        # Multiple paths founds. Plotting the graph and asking the user where to start
        Status <- 'Multiple'
      }
      
    }
    
    if(PathType == 'Long.Linear'){
      
      # Looking at all the possible diameters in the graph
      
      DiamLen <- igraph::diameter(graph = Net, directed = FALSE, unconnected = TRUE)
      
      Pattern <- igraph::graph.ring(n = DiamLen, directed = FALSE, mutual = FALSE, circular = FALSE)
      
      PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = Net, graph2 = Pattern)
      
      if (length(PossiblePaths) == 0) {
        Status <- 'None'
      }
      
      if (length(PossiblePaths) == 1) {
        Status <- 'Single'
      }
      
      if (length(PossiblePaths)> 1) {
        Status <- 'Multiple'
      }
      
      if(Circular){
        print("Circular option incompatible with Long.Linear")
      }
      
    }
    
  }
  
  
  
  if(Status == 'None'){
    print("No path found")
    return(NULL)
  }
  
  
  if(Status == 'Single'){
    print("A single path was found")
    PathToUse <- PossiblePaths[[1]]$name
  }
  
  
  if(Status == 'Multiple'){
    print(paste(length(PossiblePaths), "paths found"))
    
    TaxonList <- getTaxonMap(Results, TransfData, UseR = TRUE)
    
    InfoData <- plotPieNet(Results = PrinGraph, Graph = Net,
                           Data = TransfData, Categories = CellClass,
                           TaxonList = TaxonList,
                           NodeSizeMult = 7, ColCat = NULL, LayOut = 'circle',
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
        print("Path found")
        DONE <- TRUE
        PathToUse <- PossiblePaths[[PathToUseID]]$name
      }
      
      if(!DONE){
        print("Path not found. Please reselect the vertices")
      }
    }
  }
  
  
  # All preprocessing is done. Now we can look at gene expression over pseudotime
  
  if(Circular){
    NumericPath <- as.numeric(unlist(lapply(strsplit(rev(PathToUse), "V_"), "[[", 2)))
    NumericPath <- c(NumericPath, NumericPath[1])
  }
  
  print("Projecting cells on path")
  
  PathProjection <- OrderOnPath(PrinGraph = PrinGraph, Path = NumericPath, PointProjections = Projections)
  
  # PlotOnPath(PathProjection, StageVect)
  
  ProjectedPoints <- which(!is.na(PathProjection$PositionOnPath))
  
  SortedProj <- sort(PathProjection$PositionOnPath[ProjectedPoints], index.return=TRUE)
  
  MatGenesToPlot <- ExpressionData[ProjectedPoints[SortedProj$ix],FoundGenes]

  print("Preparing plot")
  
  if(!Circular){
    
    ggMat <- cbind(rep(SortedProj$x/sum(PathProjection$PathLen), ncol(MatGenesToPlot)),
                   as.vector(MatGenesToPlot),
                   rep(as.character(CellClass[ProjectedPoints[SortedProj$ix]]),
                       each = ncol(MatGenesToPlot)),
                   rep(colnames(MatGenesToPlot), each = nrow(MatGenesToPlot))
    )
    
    colnames(ggMat) <- c("Pseudo.Time", "Log.Gene.Exp", "Class", "Gene")
    
    ggMat <- data.frame(ggMat)
    ggMat$Pseudo.Time <- as.numeric(as.character(ggMat$Pseudo.Time))
    ggMat$Log.Gene.Exp <- as.numeric(as.character(ggMat$Log.Gene.Exp))
    
    # LM <- lm(ggMat$Log.Gene.Exp ~ ggMat$Pseudo.Time)
    # summary(LM)
    
    if(length(unique(ggMat$Class))>1){
      p <- ggplot2::ggplot(ggMat, ggplot2::aes(x = Pseudo.Time, y = Log.Gene.Exp)) +
        ggplot2::geom_point(ggplot2::aes(color = Class)) + ggplot2::facet_wrap(~ Gene, scales="free_y") +
        ggplot2::geom_smooth(color = "black", method = "loess") + ggplot2::labs(x = "Pseudo time", y = "Gene expression")
    } else {
      p <- ggplot2::ggplot(ggMat, ggplot2::aes(x = Pseudo.Time, y = Log.Gene.Exp)) +
        ggplot2::geom_point() + ggplot2::facet_wrap(~ Gene, scales="free_y") +
        ggplot2::geom_smooth(color = "black", method = "loess") + ggplot2::labs(x = "Pseudo time", y = "Gene expression")
    }
    
  } else {
    
    # We are going to include extra "shadow" points
   
    print("including virtual points")
    
    ggMat <- cbind(rep(SortedProj$x/sum(PathProjection$PathLen), ncol(MatGenesToPlot)),
                   as.vector(MatGenesToPlot),
                   rep(as.character(CellClass[ProjectedPoints[SortedProj$ix]]),
                       each = ncol(MatGenesToPlot)),
                   rep(colnames(MatGenesToPlot), each = nrow(MatGenesToPlot))
    )
    
    ggMat <- cbind(ggMat, rep("Real", nrow(ggMat)))
     
    ggMat_Minus <- ggMat[ggMat[, 1]>.9, ]
    ggMat_Minus[,5] <- "Virtual"
    ggMat_Minus[,1] <- as.numeric(ggMat_Minus[,1]) - 1
    
    ggMat_Plus <- ggMat[ggMat[, 1]<.1, ]
    ggMat_Plus[,5] <- "Virtual"
    ggMat_Plus[,1] <- as.numeric(ggMat_Plus[,1]) + 1
    
    ggMat <- rbind(ggMat, ggMat_Plus, ggMat_Minus)
    
    colnames(ggMat) <- c("Pseudo.Time", "Log.Gene.Exp", "Class", "Gene", "Status")
    
    ggMat <- data.frame(ggMat)
    ggMat$Pseudo.Time <- as.numeric(as.character(ggMat$Pseudo.Time))
    ggMat$Log.Gene.Exp <- as.numeric(as.character(ggMat$Log.Gene.Exp))
    
    # LM <- lm(ggMat$Log.Gene.Exp ~ ggMat$Pseudo.Time)
    # summary(LM)
    
    if(length(unique(ggMat$Class))>1){
      p <- ggplot2::ggplot(ggMat, ggplot2::aes(x = Pseudo.Time, y = Log.Gene.Exp)) +
        ggplot2::geom_point(ggplot2::aes(color = Class, alpha = Status)) + ggplot2::facet_wrap(~ Gene, scales="free_y") +
        ggplot2::geom_smooth(color = "black", method = "loess") + ggplot2::labs(x = "Pseudo time", y = "Gene expression") +
        ggplot2::scale_alpha_discrete(range = c(1, 0.2)) +
        ggplot2::geom_vline(xintercept = c(0,1), linetype = "dashed")
    } else {
      p <- ggplot2::ggplot(ggMat, ggplot2::aes(x = Pseudo.Time, y = Log.Gene.Exp)) +
        ggplot2::geom_point(ggplot2::aes(alpha = Status)) + ggplot2::facet_wrap(~ Gene, scales="free_y") +
        ggplot2::geom_smooth(color = "black", method = "loess") + ggplot2::labs(x = "Pseudo time", y = "Gene expression") +
        ggplot2::scale_alpha_discrete(range = c(1, 0.2)) +
        ggplot2::geom_vline(xintercept = c(0,1), linetype = "dashed")
    }
    
  }
  
  print("Plotting")
  
  ggplot2::print(p) 
  
  if (ncol(ggMat) == 5) {
    ggMat <- ggMat[ggMat[,5] == "Real", ]
    ggMat <- ggMat[, -5]
  }
  
  Smoothed <- list()
  
  for (gID in unique(ggMat[,4])) {
    Smoothed[[as.character(gID)]] <- lowess(ggMat$Log.Gene.Exp ~ ggMat$Pseudo.Time)
  }
  
  return(Smoothed)
  
  
}