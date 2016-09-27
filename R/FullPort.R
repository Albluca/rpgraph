# Configuration Functions --------------------------------------------

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

getTaxonMap <- function(Graph, Data){
  
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




# Plotting Functions --------------------------------------------

plotMSDEnergyPlot <- function(PrintGraph, Main = '', Cex.Main = .7){

  attach(PrintGraph$Report)

  plot(STEP, MSEP, type = 'b', col='blue', ylim = c(0, max(c(MSEP, ENERGY))),
       xlab = "Number of steps", ylab = "MSD, Energy", lwd = 2, pch = 1,
       main = Main, cex.main = Cex.Main)
  points(STEP, ENERGY, type = 'b', col='red', lwd = 2, pch = 2)

  legend(x = "topright",legend = c('MSD','El.Energy'),
         col = c('blue', 'red'), pch = c(1,2), lwd = 2);

  detach(PrintGraph$Report)

}

accuracyComplexityPlot <- function(PrintGraph, Main = '', Cex.Main = .7,
                                   Cex.Text = .8, Mode = 'LocMin'){

  attach(PrintGraph$Report)

  plot(FVEP, URN2, type = 'b', col='green',
       xlab = "Fraction of Explained Variance",
       ylab = "Geometrical Complexity", lwd = 2,
       pch = 1, log = 'x', main = Main, cex.main = Cex.Main)

  if(Mode == 'LocMin'){
    for(i in 2:(length(URN2)-1)){
      xp = URN2[i-1]
      x = URN2[i]
      xn = URN2[i+1]
      if(x < min(c(xp,xn))){
        diff = abs(x-(xp+xn)/2);
        if(diff>0.01){
          text(x = FVEP[i], y = URN2[i], labels = PrintGraph$Report$BARCODE[i+1],
               pos = 1, cex = Cex.Text, col = 'red')
        }
      }
    }
  }
  
  
  if(is.numeric(Mode)){
    Mode = round(Mode)
    
    text(x = FVEP[2], y = URN2[2], labels = PrintGraph$Report$BARCODE[3],
         pos = 1, cex = Cex.Text, col = 'red')
    
    text(x = FVEP[length(URN2)-1], y = URN2[length(URN2)-1],
         labels = PrintGraph$Report$BARCODE[length(URN2)],
         pos = 1, cex = Cex.Text, col = 'red')
    
    if(Mode > 2){
      
      Mode <- Mode - 1
      Step <- (length(URN2) - 2)/Mode
      
      for (i in seq(from=2+Step, to = length(URN2)-1, by = Step)) {
        text(x = FVEP[round(i)], y = URN2[round(i)],
             labels = PrintGraph$Report$BARCODE[round(i)],
             pos = 1, cex = Cex.Text, col = 'red')
      }
      
    }
    
    
  }
  
  
  detach(PrintGraph$Report)

}

plotData2D <- function(Data, PrintGraph, ScaleFunction = sqrt, NodeSizeMult=1, Col=NULL,
                       CirCol="black", LineCol="black", IdCol="blue", Main = '', Cex.Main = .7,
                       Xlab = "PC1", Ylab = "PC2"){

  
  if(is.null(Col)){
    
    Col <- rep("black", ncol(data))
    
  }

  plot(Data[,1], Data[,2], col = Col, main = Main, cex.main = Cex.Main, xlab = Xlab, ylab = Ylab)

  points(PrintGraph$Nodes[,1:2], col = CirCol, cex = NodeSizeMult*do.call(what = ScaleFunction, list(PrintGraph$NodeSize)))

  if(min(PrintGraph$Edges)==0){
    PrintGraph$Edges <- PrintGraph$Edges + 1
  }

  for(i in 1:nrow(PrintGraph$Edges)){

    arrows(x0 = PrintGraph$Nodes[PrintGraph$Edges[i,1],1], y0 = PrintGraph$Nodes[PrintGraph$Edges[i,1],2],
           x1 = PrintGraph$Nodes[PrintGraph$Edges[i,2],1], y1 = PrintGraph$Nodes[PrintGraph$Edges[i,2],2],
           length = 0, col = LineCol)

  }

  text(PrintGraph$Nodes[,1:2], labels = 1:nrow(PrintGraph$Nodes), col = IdCol, pos = 2)


}

computeMetroMapLayout <- function(PrintGraph){
  
  cat("Configuring engine .")
  
  Mml <- .jnew(class = "vdaoengine/analysis/grammars/MetroMapLayout")
  Tree <- .jnew(class = "vdaoengine/analysis/grammars/Tree")
  
  Mml$principalTree <- Tree
  
  NormEdges <- PrintGraph$Edges
  
  if(!any(NormEdges == 0)){
    NormEdges <- NormEdges - 1
  } 

  Utils <- .jnew(class = "vdaoengine/analysis/grammars/Utils")
  
  Utils$setNodesForGraph(Mml$principalTree, .jarray(.jfloat(data.matrix(PrintGraph$Nodes)),dispatch=T))
  Utils$setEdgesForGraph(Mml$principalTree, .jarray(data.matrix(NormEdges),dispatch=T))
  
  Mml$principalTree$defineStarsFromPrimitiveGraphStructure()
  
  Mml$computeLayout();
  
  NodePositions2D = .jevalArray(Mml$computedLayout$getNodePositions(), simplify = TRUE)
  
  return(NodePositions2D)
  
}

plotMappedData2D <- function(data, PrintGraph){
  
  plot(data, cex = PrintGraph$NodeSize/4)
  
  if(min(PrintGraph$Edges)==0){
    PrintGraph$Edges <- PrintGraph$Edges + 1
  }
  
  for(i in 1:nrow(PrintGraph$Edges)){
    
    arrows(x0 = data[PrintGraph$Edges[i,1],1], y0 = data[PrintGraph$Edges[i,1],2],
           x1 = data[PrintGraph$Edges[i,2],1], y1 = data[PrintGraph$Edges[i,2],2],
           length = 0)
    
  }
  
  text(data[,1:2], labels = 1:nrow(data), col = "blue", pos = 2)
  
}


plotPieNet <- function(Results, Data, Categories, LayOut = 'nicely', Main="",
                       ScaleFunction = sqrt, NodeSizeMult = 1, ColCat = NULL,
                       DirectionMat = NULL, Thr = 0.05, Arrow.size = .5, BiggestComponents = FALSE) {
  
  CatNames <- Categories
  
  require('igraph')
  require('rJava')
  
  if(is.null(ColCat)){
    ColCat <- c(rainbow(length(unique(Categories))), NA)
  } else {
    ColCat <- c(ColCat, NA)
  }
  
  TaxonList <- getTaxonMap(Graph = makeGraph(Results), Data = Data)
  TaxonCatNoNone <- list()
  
  for (i in 1:length(TaxonList)) {
    TTab <- table(CatNames[TaxonList[[i]]])
    TaxonCatNoNone[[as.character(paste("V_",i, sep = ''))]] <- TTab 
  }
  
  levels(CatNames) <- c(levels(CatNames), "None")
  
  TaxonCat <- list()
  
  for (i in 1:length(TaxonList)) {
    TTab <- table(CatNames[TaxonList[[i]]])
    if(sum(TTab) == 0){
      TTab["None"] <- 1
    }
    TaxonCat[[as.character(paste("V_",i, sep = ''))]] <- TTab 
  }
  
  
  # Generate An Igraph net
  
  library(igraph)
  
  if(min(Results$Edges) == 0){
    Results$Edges <- Results$Edges + 1
  }
  
  if(is.null(DirectionMat)){
    
    Net <- graph.empty(n = length(unique(as.vector(Results$Edges))), directed = FALSE)
    V(Net)$name <- paste("V_", unique(as.vector(Results$Edges)), sep = '')
    
    for (i in 1:nrow(Results$Edges)) {
      Net <- add.edges(graph = Net, paste("V_", Results$Edges[i,], sep = ''))
    }
    
  } else {
    
    Net <- graph.empty(n = length(unique(as.vector(Results$Edges))), directed = TRUE)
    V(Net)$name <- paste("V_", unique(as.vector(Results$Edges)), sep = '')
    
    for (i in 1:nrow(Results$Edges)) {
      
      Vertices <- paste("V_", Results$Edges[i,], sep = '')
      
      Dir1 <- which(DirectionMat$Source == Vertices[1] & DirectionMat$Target == Vertices[2])
      
      Dir2 <- which(DirectionMat$Source == Vertices[2] & DirectionMat$Target == Vertices[1])
      
      if(length(Dir1) == 0 & length(Dir2) == 0){
        print("Error: Directionality matrix incompatible")
        return(NULL)
      }
      
      if(length(Dir1) == 1){
        if(is.na(DirectionMat$P.val[Dir1])){
          next()
        }
        if(DirectionMat$Direction[Dir1] == 0 | as.numeric(DirectionMat$P.val[Dir1]) > Thr){
          Net <- add.edges(graph = Net, c(Vertices, rev(Vertices)))
          next()
        }
        if(DirectionMat$Direction[Dir1] == 1){
          Net <- add.edges(graph = Net, Vertices)
          next()
        }
        if(DirectionMat$Direction[Dir1] == 2){
          Net <- add.edges(graph = Net, rev(Vertices))
          next()
        }
      }
      
      
      if(length(Dir2) == 1){
        if(is.na(DirectionMat$P.val[Dir2])){
          next()
        }
        if(DirectionMat$Direction[Dir2] == 0 | as.numeric(DirectionMat$P.val[Dir2]) > Thr){
          Net <- add.edges(graph = Net, c(Vertices, rev(Vertices)))
          next()
        }
        if(DirectionMat$Direction[Dir2] == 1){
          Net <- add.edges(graph = Net, rev(Vertices))
          next()
        }
        if(DirectionMat$Direction[Dir2] == 2){
          Net <- add.edges(graph = Net, Vertices)
          next()
        }
        
      }
      
    }
    
  }
  
  
  OrgNetSize <- vcount(Net)
  
  if(BiggestComponents){
    
    Comp <- components(graph = Net, mode = "weak")
    BigComp <- which.max(Comp$csize) 
    
    Net <- induced.subgraph(graph = Net, vids = names(Comp$membership[Comp$membership == BigComp]))
    
  }
  
  
  PieColList <- list()
  
  for(i in 1: length(TaxonCat)){
    PieColList[[i]] <- ColCat
  }
  
  
  LayOutDONE <- FALSE
  
  if(LayOut == 'metro'){
    RestrNodes <- computeMetroMapLayout(Results)
    LayOutDONE <- TRUE
  }
  
  if(LayOut == 'tree'){
    RestrNodes <- layout_as_tree(graph = as.undirected(Net, mode = 'collapse'))
    LayOutDONE <- TRUE
  }
  
  if(LayOut == 'circle'){
    IsoGaph <- graph.ring(n = vcount(Net), directed = FALSE, circular = TRUE)
    Iso <- graph.get.isomorphisms.vf2(as.undirected(Net, mode = 'collapse'), IsoGaph)
    if(length(Iso)>0){
      VerOrder <- V(Net)[Iso[[1]]]
      RestrNodes <- layout_in_circle(graph = Net, order = VerOrder)
      LayOutDONE <- TRUE
    } else {
      LayOut = 'nicely'
    }
  }
  
  if(LayOut == 'circle_line'){
    IsoGaph <- graph.ring(n = vcount(Net), directed = FALSE, circular = FALSE)
    Iso <- graph.get.isomorphisms.vf2(as.undirected(Net, mode = 'collapse'), IsoGaph)
    if(length(Iso) > 0){
      VerOrder <- V(Net)[Iso[[1]]]
      RestrNodes <- layout_in_circle(graph = Net, order = VerOrder)
      LayOutDONE <- TRUE
    } else {
      LayOut = 'nicely'
    }
    
  }
  
  if(LayOut == 'nicely'){
    RestrNodes <- layout_nicely(graph = Net)
    LayOutDONE <- TRUE
  }
  
  if(!LayOutDONE){
    print(paste("LayOut =", LayOut, "unrecognised"))
    return(NULL)
  }
  
  Indexes <- as.integer(factor(V(Net)$name, levels = paste("V_", 1:OrgNetSize, sep='')))
  
  plot(Net, layout = RestrNodes[,1:2], main = Main,
       vertex.shape="pie", vertex.pie.color = PieColList[Indexes],
       vertex.pie=TaxonCat[Indexes], vertex.pie.border = NA,
       vertex.size=NodeSizeMult*do.call(what = ScaleFunction, list(unlist(lapply(TaxonCatNoNone[Indexes], sum)))),
       edge.color = "black", edge.arrow.size = Arrow.size)
  
  CatReord <- TaxonCat[Indexes]
  CatReordMat <- NULL
  for (i in 1:length(CatReord)) {
    CatReordMat <- rbind(CatReordMat, CatReord[[i]])
  }
  
  InfoData <- list(CatInfo = CatReordMat, ColInfo = ColCat, Net = Net, TaxonList = TaxonList)
  
  return(InfoData)
  
}









Initialize3d <- function(){
  require(rgl)
  open3d()
}




plotData3D <- function(data, PrintGraph, ScaleFunction = sqrt, NodeSizeMult=1, Col=NULL,
                       CirCol="black", LineCol="black", IdCol="blue", Main = '', Cex.Main = .7,
                       PlotProjections = FALSE, ProjectionLines = NULL,
                       Xlab = "PC1", Ylab = "PC2", Zlab = "PC3", DirectionMat = NULL, Thr = 0.05, Plot.ly = FALSE){

  # Initialize3d()
  
  if(is.null(Col)){
    Col <- rep("black", ncol(data))
  }
  
  if(is.null(Col)){
    Col <- rep("grey", ncol(data))
  }
  

  if(Plot.ly){
    
    PlotData <- TransfData[-ToRem,1:3]
    colnames(PlotData) <- c("x", "y", "z")
    PlotData <- data.frame(PlotData)
    PlotData$x <- as.numeric(as.character(PlotData$x))
    PlotData$y <- as.numeric(as.character(PlotData$y))
    PlotData$z <- as.numeric(as.character(PlotData$z))
    
    f <- list(
      family = "Courier New, monospace",
      size = 18,
      color = "#7f7f7f"
    )
    x <- list(
      title = Xlab,
      titlefont = f
    )
    y <- list(
      title = Ylab,
      titlefont = f
    )
    z <- list(
      title = Zlab,
      titlefont = f
    )
    
    
    p <- plot_ly(x = rnorm(10), y = rnorm(10), mode = "markers") %>%
      layout(xaxis = x, yaxis = y)
    p
    
    
    p <- plot_ly(data = PlotData, x = PlotData$x, y = PlotData$y, z = PlotData$z, type = "scatter3d", 
                 color = DayLabels.Factor[-ToRem], colors = unique(ColLabels), mode = "markers",
                 size = rep(1, nrow(PlotData))) %>% layout(xaxis = x, yaxis = y)
    p
    
    
    plotly_POST(p, filename = "Test3d", sharing = "public")
    
    
  } else {
    
    library(rgl)
    
    open3d()
    
    plot3d(TransfData[-ToRem,1], TransfData[-ToRem,2], TransfData[-ToRem,3], col=ColLabels[-ToRem],
           size=3, main = Main, cex.main = Cex.Main, xlab = Xlab, ylab = Ylab, zlab = Zlab, top = TRUE) 
    
    text3d(PrintGraph$Nodes[,1:3], texts = 1:nrow(data), col = IdCol)
    
    plot3d(PrintGraph$Nodes[,1:3], type = 's', radius = NodeSizeMult*do.call(what = ScaleFunction, list(PrintGraph$NodeSize)),
           add = TRUE, alpha=0.3)
    
    if(min(PrintGraph$Edges)==0){
      PrintGraph$Edges = PrintGraph$Edges + 1
    }
    
    if(is.null(DirectionMat)){
      for(i in 1:nrow(PrintGraph$Edges)){
        
        PCoords <- rbind(PrintGraph$Nodes[PrintGraph$Edges[i,1],1:3],
                         PrintGraph$Nodes[PrintGraph$Edges[i,2],1:3])
        
        plot3d(PCoords, type = 'l', add = TRUE)
        
      }
    } else {
      
      for(i in 1:nrow(DirectionMat)){
        
        SourceID <- as.integer(strsplit(DirectionMat[i, ]$Source, split = "V_")[[1]][2])
        TargetID <- as.integer(strsplit(DirectionMat[i, ]$Target, split = "V_")[[1]][2])
        Dir <- DirectionMat[i, ]$Direction
        P.val <- as.numeric(DirectionMat[i, ]$P.val)
        
        if(is.na(P.val)){
          next()
        }
        
        if(P.val > Thr | Dir == 0){
          PCoords <- rbind(PrintGraph$Nodes[SourceID,1:3],
                           PrintGraph$Nodes[TargetID,1:3])
          plot3d(PCoords, type = 'l', add = TRUE)
          next()
        }
        
        if(Dir == 1){
          arrow3d(p0 = PrintGraph$Nodes[SourceID,1:3], p1 = PrintGraph$Nodes[TargetID,1:3], type = "rotation", add=TRUE, s= .5)
          next()
        }
        
        if(Dir == 2){
          arrow3d(p1 = PrintGraph$Nodes[SourceID,1:3], p0 = PrintGraph$Nodes[TargetID,1:3], type = "rotation", add=TRUE, s= .5)
          next()
        }
        
      }
      
    }
    
    
    
    if(PlotProjections){
      
      TaxonList <- getTaxonMap(Graph = makeGraph(PrintGraph), Data = data)
      
      for(i in 1:length(TaxonList)){
        
        if(!is.na(TaxonList[[i]][1])){
          for(j in 1:length(TaxonList[[i]])){
            PCoords <- rbind(PrintGraph$Nodes[i,1:3],
                             data[TaxonList[[i]][j],1:3])
            plot3d(PCoords, type = 'l', add = TRUE, col = ProjectionLines[TaxonList[[i]][j]])
          }
        }
        
      }
      
    }
    
  }

}





# Supporting Functions --------------------------------------------


CheckArguments <- function(FunName, ArgToCheck, Filter = FALSE){

  ExtractedArgs <- do.call(args, list(name=FunName))

  AvailFunArgs <- names(as.list(ExtractedArgs))
  if(any(AvailFunArgs == "...")){
    print("Function arguments checking unavailable due to ... in its definition")
    return(ArgToCheck)
  } else {
    if(all(ArgToCheck %in% AvailFunArgs)){
      print("All arguments passed consistency check")
      return(ArgToCheck)
    } else {

      if(Filter){
        print(paste(ArgToCheck[!(ArgToCheck %in% AvailFunArgs)], "unrecognised as function argument"))
        print("They will still be passed, but please check that the execution is fine")
        return(ArgToCheck)
      } else {
        print(paste(ArgToCheck[!(ArgToCheck %in% AvailFunArgs)], "unrecognised as function argument"))
        print("They will not be passed")
        return(ArgToCheck[ArgToCheck %in% AvailFunArgs])
      }

    }

  }

}




# PCA Functions --------------------------------------------


#' Compute the PCA for a given numeric matrix usig different methodologies
#' 
#' @param DataMatrix The matrix used to compute PCA.
#' @param Components The number of principal components to be retained.
#' @param Method The method to be used to compute PCA. The current implementation accepts 'base-svd',
#' 'base-cov', 'flashPCA-eigen', 'flashPCA-svd', and 'irlba-Lanczos'
#' @return A list composed five elements rotation
#' @return Comp the rotation matrix
#' @return Center A vector containing the velue used to center the original matrix
#' 
#' 
#' @examples
#' add(1, 1)
#' add(10, 1)

SelectComputePCA <- function(DataMatrix, Components = NULL, Method = 'base-svd', ...){
  
  require("tictoc")
  
  tic()
  
  if(is.null(Components)){
    Components <- min(c(nrow(DataMatrix), ncol(DataMatrix)))/10
    Components <- ceiling(Components)
  }
  
  ExtraArgs <- list(...)

  if(ncol(DataMatrix) > nrow(DataMatrix)){

    print("The number of columns is largr than the number of rows.")
    print("This is incompatible with some PCA functions")

  }

  if(!is.numeric(Components)){

    print(paste("Invalid number of components, using", ncol(DataMatrix)))
    Components = ncol(DataMatrix)

  } else {

    if((Components > ncol(DataMatrix)) | (Components < 1)){

      print(paste("Invalid number of components, using", ncol(DataMatrix)))
      Components = ncol(DataMatrix)

    }
  }

  if(Method == 'base-svd'){
    # Using svd decomposition
    print("Using base R svd PCA")
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(prcomp, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(x = DataMatrix)
    }

    if(!any(names(FunArg)=="center")){
      FunArg <- append(FunArg, list(center=TRUE))
    }

    PCA <- do.call(what = prcomp, args = FunArg)

    RetVal <- list()
    RetVal$centers <- PCA$center
    RetVal$Comp <- PCA$rotation[,1:Components]
    RetVal$ExpVar <- (PCA$sdev[1:Components]^2)/sum(PCA$sdev^2)

    toc()
    return(RetVal)
  }


  if(Method == 'base-cov'){
    # Using svd decomposition
    print("Using base R covariance PCA")
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(princomp, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(x = DataMatrix)
    }

    if(any(names(FunArg)=="scores")){
      FunArg[["scores"]] <- TRUE
      print("The scores argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(scores=TRUE))
    }

    PCA <- do.call(what = princomp, args = FunArg)

    RetVal <- list()
    RetVal$Comp <- PCA$loadings[,1:Components]
    RetVal$centers <- PCA$center
    RetVal$ExpVar <- (PCA$sdev[1:Components]^2)/sum(PCA$sdev^2)

    toc()
    return(RetVal)
  }


  if(Method == 'flashPCA-eigen'){
    # Using covariance
    print("Using flashpcaR covariance PCA")
    if(!require("flashpcaR")){
      print("Impossible to load flashpcaR package")
      return(NULL)
    }
    if(!require("bigpca")){
      print("Impossible to load bigpca package")
      return(NULL)
    }
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(X = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(X = DataMatrix)
    }

    if(any(names(FunArg)=="method")){
      FunArg[["method"]] <- "eigen"
      print("The method argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(method="eigen"))
    }

    if(any(names(FunArg)=="ndim")){
      FunArg[["ndim"]] <- Components
      print("The ndim argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(ndim=Components))
    }
    
    if(any(names(FunArg)=="do_loadings")){
      FunArg[["do_loadings"]] <- TRUE
      print("The do_loadings argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(do_loadings=TRUE))
    }

    
    
    if(!any(names(FunArg)=="stand")){
      FunArg <- append(FunArg, list(stand="center"))
    }

    if(!any(names(FunArg)=="return_scale")){
      FunArg <- append(FunArg, list(return_scale=TRUE))
    }

    PCA <- do.call(what = flashpca, args = FunArg)

    RetVal <- list()
    RetVal$Comp <- PCA$loadings
    RetVal$centers <- PCA$center

    # Estimating explained variance 
    
    RetVal$ExpVar <- estimate.eig.vpcs(eigenv = PCA$values, M=DataMatrix, estimated = FALSE)$variance.pcs

    toc()
    return(RetVal)
  }
  
  
  if(Method == 'flashPCA-svd'){
    # Using svd decomposition
    print("Using flashpcaR svd PCA")
    if(!require("flashpcaR")){
      print("Impossible to load flashpcaR package")
      return(NULL)
    }
    if(!require("bigpca")){
      print("Impossible to load bigpca package")
      return(NULL)
    }
    print("Working")
    
    # get function arguments
    
    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(X = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(X = DataMatrix)
    }
    
    if(any(names(FunArg)=="method")){
      FunArg[["method"]] <- "svd"
      print("The method argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(method="svd"))
    }
    
    if(any(names(FunArg)=="ndim")){
      FunArg[["ndim"]] <- Components
      print("The ndim argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(ndim=Components))
    }
    
    if(any(names(FunArg)=="do_loadings")){
      FunArg[["do_loadings"]] <- TRUE
      print("The do_loadings argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(do_loadings=TRUE))
    }
    
    
    
    if(!any(names(FunArg)=="stand")){
      FunArg <- append(FunArg, list(stand="center"))
    }
    
    if(!any(names(FunArg)=="return_scale")){
      FunArg <- append(FunArg, list(return_scale=TRUE))
    }
    
    PCA <- do.call(what = flashpca, args = FunArg)
    
    RetVal <- list()
    RetVal$Comp <- PCA$loadings
    RetVal$centers <- PCA$center
    
    RetVal$ExpVar <- estimate.eig.vpcs(eigenv = PCA$values, M=DataMatrix, estimated = FALSE)$variance.pcs
    
    toc()
    return(RetVal)
  }
  
  
  if(Method == 'irlba-Lanczos'){
    # Using svd decomposition
    print("Using irlba PCA")
    if(!require("irlba")){
      print("Impossible to load irlba package")
      return(NULL)
    }
    if(!require("bigpca")){
      print("Impossible to load bigpca package")
      return(NULL)
    }
    print("Working")
    
    # get function arguments
    
    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(x = DataMatrix)
    }
    
    if(any(names(FunArg)=="n")){
      FunArg[["n"]] <- Components
      print("The n argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(n=Components))
    }
    
    
    if(!any(names(FunArg)=="center")){
      FunArg <- append(FunArg, list(center=TRUE))
    }
    
    PCA <- do.call(what = prcomp_irlba, args = FunArg)
    
    RetVal <- list()
    RetVal$Comp <- PCA$rotation
    RetVal$centers <- PCA$center
    
    # Estimating explained variance 
    
    RetVal$ExpVar <- estimate.eig.vpcs(eigenv = PCA$sdev^2, M=DataMatrix, estimated = FALSE)$variance.pcs
    
    toc()
    return(RetVal)
  }
  
  
  if(Method == 'nsprcomp'){
    # Using svd decomposition
    print("Using constrained PCA")
    if(!require("nsprcomp")){
      print("Impossible to load nsprcomp package")
      return(NULL)
    }
    if(!require("bigpca")){
      print("Impossible to load bigpca package")
      return(NULL)
    }
    print("Working")
    
    # get function arguments
    
    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(nsprcomp, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(x = DataMatrix)
    }
    
    if(any(names(FunArg)=="ncomp")){
      FunArg[["ncomp"]] <- Components
      print("The ncomp argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(ncomp=Components))
    }
    
    
    if(!any(names(FunArg)=="center")){
      FunArg <- append(FunArg, list(center=TRUE))
    }
    
    PCA <- do.call(what = nsprcomp, args = FunArg)
    
    RetVal <- list()
    RetVal$Comp <- PCA$rotation
    RetVal$centers <- PCA$center
    
    # Estimating explained variance 
    
    RetVal$ExpVar <- estimate.eig.vpcs(eigenv = PCA$sdev^2, M=DataMatrix, estimated = FALSE)$variance.pcs
    
    toc()
    return(RetVal)
  }
  
  print("Unknown methods. Please check.")
  return(NULL)
  
}





















# computeElPT --------------------------------------------





computeElPT <- function(data, NumNodes, Parameters, ...){

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

  cpg$setDataSetAsMassif(.jarray(.jfloat(data.matrix(data)),dispatch=T))
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


computeElasticPrincipalGraph <- function(data, NumNodes, Method = NULL,
                                         ReduceMethod = 'none', ReduceSize = 3, ...) {

  require("rJava")
  
  ExtraArgs <- list(...)
  
  # Center data

  mv <- apply(data, 2, mean)
  data_centered <- t(t(data) -  mv)

  # Reduce dimensions?

  RedDone <- FALSE

  if(ReduceMethod != 'none'){

    print(paste("Performing dimensionality reduction using", ReduceMethod))
    
    if(ReduceMethod == 'irlba-Lanczos' & ReduceSize == ncol(DataMatrix)){
      
      print(paste("Performing dimensionality reduction using base-svd instead of irlba-Lanczos (RTM)"))
      
    } else {
      
      print(paste("Performing dimensionality reduction using", ReduceMethod))
      
    }
    
    FunArg <- append(list(DataMatrix = data_centered, Components = ReduceSize, Method = ReduceMethod), ExtraArgs[OkArgs])
    RedData <- do.call(SelectComputePCA, FunArg)
    
    NewData <- t(t(data_centered)-RetVal$centers)%*%RetVal$Comp
    
    
    if((ReduceMethod != 'base-svd' & ReduceMethod != 'base-cov') & ReduceSize < ncol(DataMatrix)){
      print(paste(signif(100*sum(RetVal$ExpVar), 4), "% (estimated) of the original variance vas retained", sep=''))
    } else {
      print(paste(signif(100*sum(RetVal$ExpVar), 4), "% of the original variance vas retained", sep=''))
    }

  } else {
    
    NewData <- data_centered
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
                          append(list(data = NewData, NumNodes = NumNodes, Parameters = Parameters),
                                 ExtraArgs))
    return(PrintGraph)

  }

}
















# Put everything together --------------------------------------------

PlotAll <- function(PCAData, CleanData, Adjust = 0, NumNodes, Method, Lab, LayOut, Do3d = FALSE){
  
  TransfData <- as.matrix(CleanData) %*% PCAData$Comp
  
  Results <- computeElasticPrincipalGraph(data = TransfData, NumNodes = NumNodes,
                                          Method = Method)
  
  # par(mfcol=c(2,2))
  plotMSDEnergyPlot(Results)
  accuracyComplexityPlot(Results)
  
  if(Adjust>0){
    
    DistMat <- as.matrix(dist(TransfData[,1:2]))
    diag(DistMat) <- Inf
    
    ToRem <- sort(apply(DistMat, 2, min), index.return=TRUE, decreasing = TRUE)$ix[1:Adjust]
    
    plotData2D(TransfData[-ToRem,], Results, Col = ColLabels[-ToRem], Main = Lab,
               Xlab = paste("PC1 (", signif(100*PCAData$ExpVar[1], 4), "%)", sep=''),
               Ylab = paste("PC2 (", signif(100*PCAData$ExpVar[2], 4), "%)", sep=''))
  } else {
    plotData2D(TransfData, Results, Col = ColLabels, Main = Lab,
               Xlab = paste("PC1 (", signif(100*PCAData$ExpVar[1], 4), "%)", sep=''),
               Ylab = paste("PC2 (", signif(100*PCAData$ExpVar[2], 4), "%)", sep=''))
  }
  
  Cols <- plotPieNet(Results = Results, Data = TransfData, Categories = DayLabels.Factor,
                     NodeSizeMult = 3, ColCat = unique(ColLabels), LayOut = LayOut)
  
  if(LayOut == 'tree'){
    legend(x = "bottom", fill=unique(ColLabels), legend = unique(DayLabels))
  } else {
    legend(x = "center", fill=unique(ColLabels), legend = unique(DayLabels))
  }
  
  if(Do3d){
    
    if(Adjust>0){
      lotData3D(TransfData, Results, Col = ColLabels, Main = Lab,
                Xlab = paste("PC1 (", signif(100*PCAData$ExpVar[1], 4), "%)", sep=''),
                Ylab = paste("PC2 (", signif(100*PCAData$ExpVar[2], 4), "%)", sep=''),
                Zlab = paste("PC3 (", signif(100*PCAData$ExpVar[3], 4), "%)", sep=''))
    } else {
      plotData3D(TransfData[-ToRem,], Results, Col = ColLabels[-ToRem], Main = Lab,
                 Xlab = paste("PC1 (", signif(100*PCAData$ExpVar[1], 4), "%)", sep=''),
                 Ylab = paste("PC2 (", signif(100*PCAData$ExpVar[2], 4), "%)", sep=''),
                 Zlab = paste("PC3 (", signif(100*PCAData$ExpVar[3], 4), "%)", sep=''))
    }
    
  }
  
  return(Results)
  
}





# Check directionality --------------------------------------------


CheckDirectionality <- function(Data, Results, OrderClass, Depth = NULL){
  
  TaxonList <- getTaxonMap(Graph = makeGraph(Results), Data = Data)
  
  if(min(Results$Edges) == 0){
    Results$Edges <- Results$Edges + 1
  }
  
  library(igraph)
  
  Net <- graph.empty(n = length(unique(as.vector(Results$Edges))), directed = FALSE)
  V(Net)$name <- paste("V_", unique(as.vector(Results$Edges)), sep = '')
  
  for (i in 1:nrow(Results$Edges)) {
    Net <- add.edges(graph = Net, paste("V_", Results$Edges[i,], sep = ''))
  }
  
  EdgDir <- NULL
  
  if(is.null(Depth)){
    Depth <- vcount(Net)
  }
  
  for (Edg in E(Net)) {
    
    tNet <- delete_edges(Net, Edg)
    
    Verts <- V(Net)$name[get.edges(Net, Edg)]
    
    Nei <- neighborhood(graph = tNet, order = Depth, nodes = Verts)
    VertSeq1 <- unlist(strsplit(V(tNet)$name[Nei[[1]]], split = "V_"))
    VertSeq1 <- as.numeric(VertSeq1[-seq(from = 1, by = 2, along.with = VertSeq1)])
    
    VertSeq2 <- unlist(strsplit(V(tNet)$name[Nei[[2]]], split = "V_"))
    VertSeq2 <- as.numeric(VertSeq2[-seq(from = 1, by = 2, along.with = VertSeq2)])
    
    Clean_Levels <- unlist(TaxonList[VertSeq1])
    Order1 <- OrderClass[Clean_Levels[!is.na(Clean_Levels)]]
    
    Clean_Levels <- unlist(TaxonList[VertSeq2])
    Order2 <- OrderClass[Clean_Levels[!is.na(Clean_Levels)]]
    
    if(sum(is.finite(Order1)) < 3 | sum(is.finite(Order2)) < 3 ){
      EdgDir <- rbind(EdgDir, c(Verts, NA, NA))
      next()
    }
    
    WT <- wilcox.test(x = Order1, y = Order2, na.rm=TRUE)
    
    if(median(Order1, na.rm = TRUE) < median(Order2, na.rm = TRUE)){
      EdgDir <- rbind(EdgDir, c(Verts, 1, WT$p.value))
      next()
    }
    
    if(median(Order1, na.rm = TRUE) > median(Order2, na.rm = TRUE)){
      EdgDir <- rbind(EdgDir, c(Verts, 2, WT$p.value))
      next()
    }
    
    if(median(Order1, na.rm = TRUE) == median(Order2, na.rm = TRUE)){
      EdgDir <- rbind(EdgDir, c(Verts, 0, WT$p.value))
      next()
    }
    
  }
  
  colnames(EdgDir) <- c("Source", "Target", "Direction", "P.val")
  
  return(data.frame(EdgDir, stringsAsFactors = FALSE))
  
}


GetDirectionalityIndex <- function(InfoData) {
  
  VertDist <- distances(InfoData$Net, mode = 'out')
  Dists <- cbind(which(lower.tri(VertDist), arr.ind = TRUE), VertDist[lower.tri(VertDist)], t(VertDist)[lower.tri(t(VertDist))])
  Dists <- cbind(V(InfoData$Net)$name[Dists[,1]], V(InfoData$Net)$name[Dists[,2]], Dists[,3:4])
  colnames(Dists) <- c("V1", "V2", "V1->V2", "V2->V1")
  
  DirectionalityIndex <- (sum(is.infinite(as.numeric(Dists[,3])) &
                                !is.infinite(as.numeric(Dists[,4]))) + 
                            sum(!is.infinite(as.numeric(Dists[,3])) &
                                  is.infinite(as.numeric(Dists[,4]))))/sum(!is.infinite(as.numeric(Dists[,3])) |
                                                                             !is.infinite(as.numeric(Dists[,4])))
  
  return(DirectionalityIndex)
}






# CompareNets -------------------------------------------------------------


CompareNet <- function(G1, G2, RemNodes = 2, Tries = 10000, DoIso = FALSE) {
  
  if(DoIso){
    Full_Iso <- graph.get.isomorphisms.vf2(graph1 = G1, graph2 = G2)
    
    if(length(Full_Iso)>0){
      return(0)
    }
  }
  
  
  pb <- txtProgressBar(min = 1, max = Tries, initial = 1, style = 3)
  
  for (Retries in 1:Tries) {
    setTxtProgressBar(pb, Retries)
    RemVert <- sample(x = V(G1), size = RemNodes, replace = FALSE)
    tNet <- delete_vertices(G1, RemVert)
    Part_Iso <- graph.get.subisomorphisms.vf2(graph1 = G2, graph2 = tNet)
    if(length(Part_Iso)>0){
      RetVal <- list(rem)
      close(pb)
      return(list(SubIso = Part_Iso, RemVert = RemVert))
    }
  }
  
  close(pb)
  return(NULL)
  
}




