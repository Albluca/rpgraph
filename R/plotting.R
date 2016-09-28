
# Plotting Functions (Diagnostic) --------------------------------------------

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
                                   Cex.Text = .8, Mode = 'LocMin', Xlims = c(0, 1)){

  attach(PrintGraph$Report)

  plot(FVEP, URN2, type = 'b', col='green',
       xlab = "Fraction of Explained Variance",
       ylab = "Geometrical Complexity", lwd = 2,
       pch = 1, main = Main, cex.main = Cex.Main,
       xlim = Xlims)

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

# Plotting Functions (2D plots) --------------------------------------------

plotData2D <- function(Data, PrintGraph, PCAData=NULL, ScaleFunction = sqrt, NodeSizeMult=1, Col=NULL,
                       CirCol="black", LineCol="black", IdCol="blue", Main = '', Cex.Main = .7,
                       Xlab = "", Ylab = ""){

  Data <- data.matrix(Data)

  if(is.null(PCAData)){
    print("Computing PCA in order to project data on 2 dimensions")
    if(ncol(Data) > 5){
      PCAData <- SelectComputePCA(DataMatrix = Data, Components = 2, Method = 'irlba-Lanczos', center=TRUE, scale.=TRUE)
    } else {
      PCAData <- SelectComputePCA(DataMatrix = Data, Components = 2, Method = 'base-svd', center=TRUE, scale.=TRUE)
    }

    # centering data
    TransfData <- Data - PCAData$centers

    # Scaling data
    TransfData <- TransfData/PCAData$scale

    # Rotating data
    TransfData <- TransfData %*% PCAData$Comp

    if(Xlab == ""){
      Xlab = paste("PC1 (", signif(100*PCAData$ExpVar[1], 4), "%)", sep='')
    }

    if(Ylab == ""){
      Ylab = paste("PC2 (", signif(100*PCAData$ExpVar[2], 4), "%)", sep='')
    }

  } else {

    # centering data
    TransfData <- Data - PCAData$centers

    # Scaling data
    TransfData <- TransfData/PCAData$scale

    # Rotating data
    TransfData <- TransfData %*% PCAData$Comp

    if(Xlab == ""){
      Xlab = paste("PC1 (", signif(100*PCAData$ExpVar[1], 4), "%)", sep='')
    }

    if(Ylab == ""){
      Ylab = paste("PC2 (", signif(100*PCAData$ExpVar[2], 4), "%)", sep='')
    }

  }

  if(is.null(Col)){
    Col <- rep("gray", nrow(Data))
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






plotPieNet <- function(Results, Data, Categories, Graph = NULL, TaxonList = NULL, LayOut = 'nicely', Main="",
                       ScaleFunction = sqrt, NodeSizeMult = 1, ColCat = NULL,
                       DirectionMat = NULL, Thr = 0.05, Arrow.size = .5, BiggestComponents = FALSE) {

  CatNames <- Categories

  if(is.null(ColCat)){
    ColCat <- c(rainbow(length(unique(Categories))), NA)
  } else {
    ColCat <- c(ColCat, NA)
  }

  if(is.null(TaxonList)){
    print("TaxonList will be computed. Consider do that separetedly")
    TaxonList <- getTaxonMap(Graph = makeGraph(Results), Data = Data)
  }

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

  if(is.null(Graph)){
    print("A graph will be constructed. Consider do that separatedly")
  } else {

    # Generate An Igraph net

    Net <- ConstructGraph(Results = Results, DirectionMat = DirectionMat, Thr = Thr)

  }


  OrgNetSize <- igraph::vcount(Net)

  if(BiggestComponents){

    Comp <- igraph::components(graph = Net, mode = "weak")
    BigComp <- which.max(Comp$csize)

    Net <- igraph::induced.subgraph(graph = Net, vids = names(Comp$membership[Comp$membership == BigComp]))

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
    RestrNodes <- igraph::layout_as_tree(graph = igraph::as.undirected(Net, mode = 'collapse'))
    LayOutDONE <- TRUE
  }

  if(LayOut == 'circle'){
    IsoGaph <- igraph::graph.ring(n = igraph::vcount(Net), directed = FALSE, circular = TRUE)
    Iso <- igraph::graph.get.isomorphisms.vf2(igraph::as.undirected(Net, mode = 'collapse'), IsoGaph)
    if(length(Iso)>0){
      VerOrder <- V(Net)[Iso[[1]]]
      RestrNodes <- igraph::layout_in_circle(graph = Net, order = VerOrder)
      LayOutDONE <- TRUE
    } else {
      LayOut = 'nicely'
    }
  }

  if(LayOut == 'circle_line'){
    IsoGaph <- igraph::graph.ring(n = igraph::vcount(Net), directed = FALSE, circular = FALSE)
    Iso <- igraph::graph.get.isomorphisms.vf2(aigraph::s.undirected(Net, mode = 'collapse'), IsoGaph)
    if(length(Iso) > 0){
      VerOrder <- V(Net)[Iso[[1]]]
      RestrNodes <- igraph::layout_in_circle(graph = Net, order = VerOrder)
      LayOutDONE <- TRUE
    } else {
      LayOut = 'nicely'
    }

  }

  if(LayOut == 'nicely'){
    RestrNodes <- igraph::layout_nicely(graph = Net)
    LayOutDONE <- TRUE
  }

  if(!LayOutDONE){
    print(paste("LayOut =", LayOut, "unrecognised"))
    return(NULL)
  }

  Indexes <- as.integer(factor(V(Net)$name, levels = paste("V_", 1:OrgNetSize, sep='')))

  igraph::plot(Net, layout = RestrNodes[,1:2], main = Main,
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






# Plotting Functions (3D plots) --------------------------------------------


Initialize3d <- function(){
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


