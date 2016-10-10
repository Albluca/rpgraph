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
                                 Path = 'ask', PathType = 'Long.Linear') {
  
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
  
  Net <- ConstructGraph(Results = PrinGraph, DirectionMat = NULL)
  
  Status <- 'unDone'
  
  if(Path != 'ask'){
    
    igraph::induced.subgraph(graph = Net, vids = PathToUse)
    
  }
  
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
    
  }
  
  
  if(Status == 'None'){
    print("No path found")
    return(NULL)
  }
  
  
  if(Status == 'Single'){
    print("A single path was found")
    PathToUseID <- PossiblePaths[[1]]$name
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
        PathToUseID <- PossiblePaths[[PathToUseID]]$name
      }
      
      if(!DONE){
        print("Path not found. Please reselect the vertices")
      }
    }
  }
  
  
  
  
  
  
  
}