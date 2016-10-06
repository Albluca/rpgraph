
# Plotting Functions (Diagnostic) --------------------------------------------

#' Title
#'
#' @param PrintGraph 
#' @param Main 
#' @param Cex.Main 
#'
#' @return
#' @export
#'
#' @examples
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







#' Title
#'
#' @param PrintGraph 
#' @param Main 
#' @param Cex.Main 
#' @param Cex.Text = .8
#' @param Mode = 'LocMin'
#' @param Xlims = NULL
#'
#' @return
#' @export
#'
#' @examples
accuracyComplexityPlot <- function(PrintGraph, Main = '', Cex.Main = .7,
                                   Cex.Text = .8, Mode = 'LocMin', Xlims = NULL){

  attach(PrintGraph$Report)
  
  if(is.null(Xlims)){
    Xlims <- range(FVEP)
  }

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


#' Title
#'
#'
#' @importFrom plotly %>%
#' @param Data 
#' @param PrintGraph 
#' @param GroupsLab 
#' @param ScaleFunction 
#' @param NodeSizeMult 
#' @param Col 
#' @param CirCol 
#' @param ColLabels 
#' @param LineCol 
#' @param IdCol 
#' @param Main 
#' @param Cex.Main 
#' @param Xlab 
#' @param Ylab 
#' @param Plot.ly 
#'
#' @return
#' @export
#'
#' @examples
plotData2D <- function(Data, PrintGraph, GroupsLab, ScaleFunction = sqrt,
                       NodeSizeMult=1, Col=NULL, CirCol="black", ColLabels = NULL,
                       LineCol="black", IdCol="blue", Main = '', Cex.Main = .7,
                       Xlab = "PC1", Ylab = "PC2", Plot.ly = FALSE){
  
  Data <- data.matrix(Data)
  
  if(is.null(Col)){
    Col <- rainbow(length(unique(GroupsLab)))[as.integer(factor(GroupsLab))]
  }
  
  if(min(PrintGraph$Edges)==0){
    PrintGraph$Edges <- PrintGraph$Edges + 1
  }
  
  if(Plot.ly){
    
    PlotData1 <- Data[,1:2]
    
    if(is.null(rownames(PlotData1))){
      rownames(PlotData1) <- paste("R_", 1:nrow(PlotData1), sep= '')
    }
    
    PlotData2 <- PrintGraph$Nodes[,1:2]
    rownames(PlotData2) <- paste("V_", 1:nrow(PlotData2), sep='')
    PlotData3 <- c(GroupsLab, rep("Graph", nrow(PlotData2)))
    
    PlotData <- cbind(rbind(PlotData1, PlotData2), PlotData3)
    
    colnames(PlotData) <- c("x", "y", "color")
    PlotData <- data.frame(PlotData)
    PlotData$x <- as.numeric(as.character(PlotData$x))
    PlotData$y <- as.numeric(as.character(PlotData$y))
    PlotData$color <- factor(PlotData$color, levels = c(unique(GroupsLab), "Graph"))
    
    
    p <- plotly::plot_ly(x = PlotData$x, y = PlotData$y,
                 type = "scatter", mode = "markers", text = rownames(PlotData),
                 color = PlotData$color,
                 colors = c(unique(Col), "black"),
                 size = 1, sizes = c(1, 10), hoverinfo = 'text')
    
    p <- p %>% plotly::layout(xaxis = list(title = Xlab), yaxis = list(title = Ylab), title=Main)

    for(i in 1:nrow(PrintGraph$Edges)){
      
      p <- p %>% plotly::add_trace(x = PrintGraph$Nodes[PrintGraph$Edges[i,1:2],1],
                           y = PrintGraph$Nodes[PrintGraph$Edges[i,1:2],2],
                           color = PlotData$color[length(PlotData$color)], text = '', size = 0.2, mode="lines",
                           showlegend = FALSE)
    }

    p
    
  } else {
    
    plot(Data[,1], Data[,2], col = Col, main = Main, cex.main = Cex.Main, xlab = Xlab, ylab = Ylab)
    
    points(PrintGraph$Nodes[,1:2], col = CirCol, cex = NodeSizeMult*do.call(what = ScaleFunction, list(PrintGraph$NodeSize)))
    
    for(i in 1:nrow(PrintGraph$Edges)){
      
      arrows(x0 = PrintGraph$Nodes[PrintGraph$Edges[i,1],1], y0 = PrintGraph$Nodes[PrintGraph$Edges[i,1],2],
             x1 = PrintGraph$Nodes[PrintGraph$Edges[i,2],1], y1 = PrintGraph$Nodes[PrintGraph$Edges[i,2],2],
             length = 0, col = LineCol)
      
    }
    
    text(PrintGraph$Nodes[,1:2], labels = 1:nrow(PrintGraph$Nodes), col = IdCol, pos = 2)
    
  }
  

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






#' Title
#'
#' @param Results 
#' @param Data 
#' @param Categories 
#' @param Graph 
#' @param TaxonList 
#' @param LayOut 
#' @param Main 
#' @param ScaleFunction 
#' @param NodeSizeMult 
#' @param ColCat 
#' @param DirectionMat 
#' @param Thr 
#' @param Arrow.size 
#' @param BiggestComponents 
#'
#' @return
#' @export
#'
#' @examples
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
    print("TaxonList will be computed. Consider doing that separetedly")
    TaxonList <- getTaxonMap(Results = Results, Data = Data)
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
    Net <- ConstructGraph(Results = Results, DirectionMat = DirectionMat, Thr = Thr)
  } else {
    Net <- Graph
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
      VerOrder <- igraph::V(Net)[Iso[[1]]]
      RestrNodes <- igraph::layout_in_circle(graph = Net, order = VerOrder)
      LayOutDONE <- TRUE
    } else {
      Net1 <- ConstructGraph(Results = Results, DirectionMat = NULL)
      IsoGaph <- igraph::graph.ring(n = igraph::vcount(Net1), directed = FALSE, circular = TRUE)
      Iso <- igraph::graph.get.isomorphisms.vf2(igraph::as.undirected(Net1, mode = 'collapse'), IsoGaph)
      VerOrder <- igraph::V(Net1)[Iso[[1]]]
      RestrNodes <- igraph::layout_in_circle(graph = Net, order = VerOrder$name)
      LayOutDONE <- TRUE
    }
  }

  if(LayOut == 'circle_line'){
    IsoGaph <- igraph::graph.ring(n = igraph::vcount(Net), directed = FALSE, circular = FALSE)
    Iso <- igraph::graph.get.isomorphisms.vf2(igraph::as.undirected(Net, mode = 'collapse'), IsoGaph)
    if(length(Iso) > 0){
      VerOrder <- igraph::V(Net)[Iso[[1]]]
      RestrNodes <- igraph::layout_in_circle(graph = Net, order = VerOrder)
      LayOutDONE <- TRUE
    } else {
      Net1 <- ConstructGraph(Results = Results, DirectionMat = NULL)
      IsoGaph <- igraph::graph.ring(n = igraph::vcount(Net1), directed = FALSE, circular = FALSE)
      Iso <- igraph::graph.get.isomorphisms.vf2(aigraph::s.undirected(Net1, mode = 'collapse'), IsoGaph)
      VerOrder <- igraph::V(Net1)[Iso[[1]]]
      RestrNodes <- igraph::layout_in_circle(graph = Net, order = VerOrder$name)
      LayOutDONE <- TRUE
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

  Indexes <- as.integer(factor(igraph::V(Net)$name, levels = paste("V_", 1:OrgNetSize, sep='')))

  igraph::plot.igraph(Net, layout = RestrNodes[,1:2], main = Main,
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



#' Plot dqta in 3 dimensions
#'
#' @importFrom plotly %>%
#' @param Data The Data to be plotted. Each column represents a dimension and each row represents a point
#' @param PrintGraph A princial graph structure obtained by \link{computeElasticPrincipalGraph}
#' @param GroupsLab A vector of labels to indicate which group each point belongs to 
#' @param ScaleFunction A scaling function to decide the size of the spheres that represent the nodes of the principal graph
#' @param NodeSizeMult A scaling vector to decide the size of the spheres that represent the nodes of the principal graph
#' @param Col A vector of colors for the points. If NULL the cplors will be computed automatically
#' @param CirCol The color of the Spheres (it will be partially transparent)
#' @param LineCol The colors of the edtges of the graph
#' @param IdCol The color of the labels of the nodes of the graph
#' @param Main The title of the plot
#' @param Cex.Main The multiplier associated with the title of the plot
#' @param PlotProjections Should the line projecting points on the nodes to be plotted?
#' @param ProjectionLines A vector of colors for the projection lines (used only if Plot.ly = FALSE)
#' @param TaxonList A list of associations between points and nodes produced by \link{getTaxonMap}
#' @param Xlab The label of the x axis
#' @param Ylab The label of the y axis
#' @param Zlab The label of the z axis
#' @param DirectionMat A directionality structure produced by \link{CheckDirectionality}
#' @param Thr A threshold to be used for directionality reconstruction
#' @param Plot.ly A boolen indicating if Plot.ly (TRUE) or rgl (FALSE) should be used
#'
#' @export
#'
#' @examples
plotData3D <- function(Data, PrintGraph, GroupsLab, ScaleFunction = sqrt, NodeSizeMult=1,
                       Col=NULL,
                       CirCol="black", LineCol="black", IdCol="blue", Main = '', Cex.Main = .7,
                       PlotProjections = FALSE, ProjectionLines = NULL, TaxonList = NULL,
                       Xlab = "PC1", Ylab = "PC2", Zlab = "PC3", DirectionMat = NULL,
                       Thr = 0.05, Plot.ly = FALSE){
  
  Data <- data.matrix(Data)
  
  if(is.null(Col)){
    Col <- rainbow(length(unique(GroupsLab)))[as.integer(factor(GroupsLab))]
  }
  
  if(min(PrintGraph$Edges)==0){
    PrintGraph$Edges = PrintGraph$Edges + 1
  }
  
  if(Plot.ly){
    
    PlotData1 <- Data[,1:3]
    
    if(is.null(rownames(PlotData1))){
      rownames(PlotData1) <- paste("R_", 1:nrow(PlotData1), sep= '')
    }
    
    PlotData2 <- PrintGraph$Nodes[,1:3]
    rownames(PlotData2) <- paste("V_", 1:nrow(PlotData2))
    PlotData3 <- c(GroupsLab, rep("Graph", nrow(PlotData2)))
    
    PlotData <- cbind(rbind(PlotData1, PlotData2), PlotData3)
    
    colnames(PlotData) <- c("x", "y", "z", "color")
    PlotData <- data.frame(PlotData)
    PlotData$x <- as.numeric(as.character(PlotData$x))
    PlotData$y <- as.numeric(as.character(PlotData$y))
    PlotData$z <- as.numeric(as.character(PlotData$z))
    
    
    p <- plotly::plot_ly(x = PlotData$x, y = PlotData$y, z = PlotData$z,
                 type = "scatter3d", mode = "markers", text = rownames(PlotData),
                 color = c(as.character(GroupsLab), rep("Graph", nrow(PlotData2))),
                 colors = c(unique(Col), col=CirCol),
                 size = 1, sizes = c(1, 10), hoverinfo = 'text') %>%
      plotly::layout(title = Main,
             scene = list(
               xaxis = list(title = Xlab), 
               yaxis = list(title = Ylab), 
               zaxis = list(title = Zlab)))
    
    
    for(i in 1:nrow(PrintGraph$Edges)){
      
      p <- p %>% plotly::add_trace(x = PrintGraph$Nodes[PrintGraph$Edges[i,1:2],1],
                           y = PrintGraph$Nodes[PrintGraph$Edges[i,1:2],2],
                           z = PrintGraph$Nodes[PrintGraph$Edges[i,1:2],3],
                           color = "Graph", text = '', size = 1,
                           sizes = c(1, 10), mode="lines",
                           showlegend = FALSE)
    }
    
    
    if(PlotProjections){
      
      if(is.null(TaxonList)){
        print("TaxonList will be computed. Consider do that separetedly")
        TaxonList <- getTaxonMap(Graph = makeGraph(PrintGraph), Data = Data)
      }
      
      for(i in 1:length(TaxonList)){
        
        if(!is.na(TaxonList[[i]][1])){
          for(j in 1:length(TaxonList[[i]])){
            
            PCoords <- rbind(PrintGraph$Nodes[i,1:3],
                             data[TaxonList[[i]][j],1:3])
            
            p <- p %>% plotly::add_trace(x = PCoords[,1],
                                 y = PCoords[,2],
                                 z = PCoords[,3],
                                 color = as.character(GroupsLab)[TaxonList[[i]][j]], text = '', size = .5,
                                 sizes = c(1, 10), mode="lines",
                                 showlegend = FALSE, opacity = 0.5)
          }
        }
        
      }
      
    }
    
    p
    
    
  } else {
    
    rgl::plot3d(Data[,1], Data[,2], Data[,3], col=Col,
           size=3, main = Main, cex.main = Cex.Main, xlab = Xlab, ylab = Ylab, zlab = Zlab, top = TRUE) 
    
    rgl::text3d(PrintGraph$Nodes[,1:3], texts = 1:nrow(Data), col = IdCol)
    
    rgl::plot3d(PrintGraph$Nodes[,1:3], type = 's', radius = NodeSizeMult*do.call(what = ScaleFunction, list(PrintGraph$NodeSize)),
           add = TRUE, alpha=0.3, col=CirCol)
    
    
    if(is.null(DirectionMat)){
      for(i in 1:nrow(PrintGraph$Edges)){
        
        PCoords <- rbind(PrintGraph$Nodes[PrintGraph$Edges[i,1],1:3],
                         PrintGraph$Nodes[PrintGraph$Edges[i,2],1:3])
        
        rgl::plot3d(PCoords, type = 'l', add = TRUE, col=LineCol)
        
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
          rgl::plot3d(PCoords, type = 'l', add = TRUE, col=LineCol)
          next()
        }
        
        if(Dir == 1){
          rgl::arrow3d(p0 = PrintGraph$Nodes[SourceID,1:3], p1 = PrintGraph$Nodes[TargetID,1:3], type = "rotation", add=TRUE, s= .5, col=LineCol)
          next()
        }
        
        if(Dir == 2){
          rgl::arrow3d(p1 = PrintGraph$Nodes[SourceID,1:3], p0 = PrintGraph$Nodes[TargetID,1:3], type = "rotation", add=TRUE, s= .5, col=LineCol)
          next()
        }
        
      }
      
    }
    
    
    
    if(PlotProjections){
      
      if(is.null(TaxonList)){
        print("TaxonList will be computed. Consider do that separetedly")
        TaxonList <- getTaxonMap(Graph = makeGraph(PrintGraph), Data = Data)
      }
      
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

