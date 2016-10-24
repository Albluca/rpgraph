# A function to standardize the Initial steps of the analysis

#' Title
#'
#' @param ExpressionMatrix 
#' @param Grouping 
#' @param GeneSet 
#' @param StageAssociation 
#' @param TopQ 
#' @param LowQ 
#' @param AggressiveStaging 
#' @param PathOpt 
#' @param GeneOpt 
#' @param GeneDetectedFilter 
#' @param GeneCountFilter 
#' @param MinCellExp 
#' @param VarFilter 
#' @param LogTranform 
#' @param Centering 
#' @param Scaling 
#' @param nDim 
#' @param nPoints 
#' @param Data.Return 
#' @param GeneToPlot 
#'
#' @return
#' @export
#'
#' @examples
StudyCellCycles <- function(ExpressionMatrix, Grouping, GeneSet = NULL,
                            StageAssociation = NULL, TopQ = .9, LowQ = .1, AggressiveStaging = FALSE,
                            PathOpt = "switch", GeneOpt = 10,
                            GeneDetectedFilter = 2.5, GeneCountFilter = 2.5,
                            MinCellExp = 1, VarFilter = 0, LogTranform = TRUE, Centering = FALSE,
                            Scaling = FALSE, nDim = NULL, nPoints = 20,
                            Data.Return = FALSE, GeneToPlot) {
  
  print(paste("Expression matrix contains", nrow(ExpressionMatrix), "cells and", nrow(ExpressionMatrix), "genes"))
  
  # Cell filtering ---------------------------------------------------------------
  
  print("Stage I - Cell filtering")
  
  hist(apply(ExpressionMatrix > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
       freq = TRUE, ylab = "Number of cells")
  
  OutExpr <- scater::isOutlier(rowSums(ExpressionMatrix>0), nmads = GeneDetectedFilter)
  
  
  hist(apply(ExpressionMatrix, 1, sum), main = "Reads per cell", xlab = "Reads count",
       freq = TRUE, ylab = "Number of cells")
  
  OutCount <- scater::isOutlier(rowSums(ExpressionMatrix>0), nmads = GeneCountFilter)

  
  print(paste(sum(OutExpr & OutCount), "Cells will be removed due to both filtering"))
  print(paste(sum(OutExpr), "Cells will be removed due to gene count filtering"))
  print(paste(sum(OutCount), "Cells will be removed due to read count filtering"))
  print(paste(sum(!(OutExpr & OutCount)), "Cells will be used for analysis")) 

  Grouping <- Grouping[!(OutExpr & OutCount)]
  NormExpressionMatrix <- ExpressionMatrix[!(OutExpr & OutCount),]
  
  readline("Press any key")
  
  par(mfcol=c(1,2))
  
  hist(apply(NormExpressionMatrix > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
       freq = TRUE, ylab = "Number of cells (After filtering)")
  
  hist(apply(NormExpressionMatrix, 1, sum), main = "Reads per cell", xlab = "Reads count",
       freq = TRUE, ylab = "Number of cells (After filtering)")
  
  
  # Gene filtering ---------------------------------------------------------------
  
  print("Stage II - Genes filtering")
  
  print(paste("Removing genes that are detected in less than", MinCellExp, "cell(s)"))
  
  GeneFilterBool <- apply(NormExpressionMatrix > 0, 2, sum) < MinCellExp
  NormExpressionMatrix <- NormExpressionMatrix[, !GeneFilterBool]

  print(paste(sum(GeneFilterBool), "genes will be removed"))
  print(paste(ncol(NormExpressionMatrix), "genes are available for the analysis"))
  
  if(LogTranform){
    print("Transforming using pseudo counts (Log10(x+1))")
    NormExpressionMatrix <- log10(NormExpressionMatrix + 1)
  }
  
  print("Plotting variance distribution (For reference only)")
    
  par(mfcol=c(1,2))

  VarVect <- apply(NormExpressionMatrix, 2, var)
  
  plot(density(VarVect, from=min(VarVect)), xlab = "Variance", main = "Variance distribution")
  abline(v=quantile(VarVect, c(.01, .05, .25, .50)), lty=2, col=c('red', 'green', 'blue', "black"))
  legend(x = "topright", legend = c("1% Q", "5% Q", "25% Q", "50% Q"),
         col=c('red', 'green', 'blue', "black"), lty=2)
  
  plot(density(VarVect, from=min(VarVect)), xlab = "Variance (1% - 75% Q)",
       xlim = quantile(VarVect, c(.01, .75)), main = "Variance distribution")
  abline(v=quantile(VarVect, c(.01, .05, .25, .50)), lty=2, col=c('red', 'green', 'blue', "black"))
    
  legend(x = "topright", legend = c("1% Q", "5% Q", "25% Q", "50% Q"),
         col=c('red', 'green', 'blue', "black"), lty=2)
  
  
  if(!is.null(GeneSet)){
    SelGenes <- intersect(colnames(NormExpressionMatrix), GeneSet)
    
    print(paste(length(SelGenes), "selected for analysis"))
    
    if(length(SelGenes)==0){
      return()
    }
    
    NormExpressionMatrix <- NormExpressionMatrix[,SelGenes]
    
  }
  
  readline("Press any key")
  
  # Gene Scaling ---------------------------------------------------------------
  
  print("Stage III - Gene scaling")
  
  if(Centering){
    print("Gene expression will be centered")
  } 
  
  if(Scaling){
    print("Gene expression will be scaled")
  } 
  
  NormExpressionMatrix <- scale(NormExpressionMatrix, Centering, Scaling)
  
  print("Stage IV - PCA")
  
  if(is.null(nDim)){
    nDim <- min(dim(NormExpressionMatrix))
  }
  
  NormExpressionMatrixPCA <- SelectComputePCA(NormExpressionMatrix,
                                              Components = nDim, Method = 'base-svd',
                                              center = FALSE, scale.=FALSE)
  
  RotatedExpression <- NormExpressionMatrix %*% NormExpressionMatrixPCA$Comp
  
  # Curve Fit ---------------------------------------------------------------
  
  print("Stage V - Curve fit")
  
  Results <- computeElasticPrincipalGraph(Data = RotatedExpression, NumNodes = nPoints, Method = 'CircleConfiguration')
  
  par(mfcol=c(1,2))
  
  accuracyComplexityPlot(Results, AdjFactor = 1, Mode = 'LocMin')
  plotMSDEnergyPlot(Results)
  
  
  ColCells <- rainbow(length(unique(Grouping)))
  Col = ColCells[as.integer(factor(ColCells))]
  
  par(mfcol=c(1,1))
  
  plotData2D(Data = RotatedExpression, PrintGraph = Results, Col = "black", NodeSizeMult = 0.1,
             Main = "All genes", Plot.ly = FALSE, GroupsLab = NULL,
             Xlab = paste("PC1 (", signif(100*NormExpressionMatrixPCA$ExpVar[1], 4), "%)", sep=''),
             Ylab = paste("PC2 (", signif(100*NormExpressionMatrixPCA$ExpVar[2], 4), "%)", sep=''))
  
  
  TaxonList <- getTaxonMap(Results, RotatedExpression, UseR = TRUE)
  
  ProjPoints <- projectPoints(Results = Results, Data = RotatedExpression, TaxonList=TaxonList,
                              UseR = TRUE,
                              method = 'PCALin', Dims = NULL)
  
  arrows(x0 = RotatedExpression[,1], y0 = RotatedExpression[,2],
         x1 = ProjPoints$PointsOnEdgesCoords[,1], y1 = ProjPoints$PointsOnEdgesCoords[,2], length = 0)
  
  Net <- ConstructGraph(Results = Results, DirectionMat = NULL)
  
  readline("Press any key")
  
  # Path Selection ---------------------------------------------------------------
  
  print("Stage V - Path Selection")
  
  # Start be selecting the first available path.
  
  Pattern <- igraph::graph.ring(n = nPoints, directed = FALSE, mutual = FALSE, circular = FALSE)
  
  PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = Net, graph2 = Pattern)
  
  UsedPath <- PossiblePaths[[1]]$name
  
  print("Using the following reference path")
  print(UsedPath)
  
  NodeOnGenes <- t(t(Results$Nodes %*% t(NormExpressionMatrixPCA$Comp)))
  
  NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
  

  NodeOnGenesOnPath <- NodeOnGenes[NumericPath,]
  
  
  if(is.list(StageAssociation)){
    
    print("Gene/Stage information found. Trying to Optimize")
    
    StageMatU <- NULL
    StageMatD <- NULL
    GeneCount <- rep(0, length(StageAssociation$Stages))
    
    print("Stage V.I - Associating peacks and valleys")
    
    for (Stage in 1:length(StageAssociation$Stages)) {
      
      if(exists(paste("S", Stage, "_U", sep = ""), where=StageAssociation)) {
        
        StageGenes <- unlist(StageAssociation[paste("S", Stage, "_U", sep = "")], use.names = FALSE)
        
        if(length(intersect(StageGenes, colnames(NodeOnGenesOnPath)))>0){
          StageTracks <- NodeOnGenesOnPath[, StageGenes]
          SignStageMat <- sign(t(StageTracks) - apply(StageTracks, 2, quantile, TopQ))
          
          SignStageMat[SignStageMat <= 0] <- NA
          SignStageMat[SignStageMat > 0] <- Stage
          
          StageMatU <- rbind(StageMatU, SignStageMat)
          GeneCount[Stage] <- GeneCount[Stage] + ncol(StageTracks)
        }
        
      }
      
      if(exists(paste("S", Stage, "_D", sep = ""), where=StageAssociation)) {
        
        StageGenes <- unlist(StageAssociation[paste("S", Stage, "_U", sep = "")], use.names = FALSE)
        
        if(length(intersect(StageGenes, colnames(NodeOnGenesOnPath)))>0){
          StageTracks <- NodeOnGenesOnPath[, StageGenes]
          SignStageMat <- sign(t(StageTracks) - apply(StageTracks, 2, quantile, LowQ))
          
          SignStageMat[SignStageMat >= 0] <- NA
          SignStageMat[SignStageMat < 0] <- Stage
          
          StageMatD <- rbind(StageMatD, SignStageMat)
          GeneCount[Stage] <- GeneCount[Stage] + ncol(StageTracks)
        }
        
      }
      
    }
    
    
    SummaryStageMat <- NULL
    
    for (Stage in 1:length(StageAssociation$Stages)) {
      
      StageCount <- rep(0, nPoints)
      
      if(!is.null(StageMatU)){
        dim(StageMatU) <- c(length(StageMatU)/nPoints, nPoints)
        StageCount <- StageCount + colSums(StageMatU == Stage, na.rm = TRUE)
      }
      
      if(!is.null(StageMatD)){
        dim(StageMatD) <- c(length(StageMatD)/nPoints, nPoints)
        StageCount <- StageCount + colSums(StageMatD == Stage, na.rm = TRUE)
      }
      
      SummaryStageMat <- rbind(SummaryStageMat, StageCount)
      
    }
    
    
    print("Stage V.II - Maximising stage association")
    
    SummaryStageMat <- SummaryStageMat/GeneCount
    
    StageNodeAssociation <- apply(SummaryStageMat, 2, which.max)
    
    NodeSize <- unlist(lapply(lapply(TaxonList, is.finite), sum))
    
    StageNodeAssociationSmooth <- StageNodeAssociation
    
    for(i in 1:length(NodeSize)){
      
      SmoothedAssociation <- SmoothFilter(CircShift(StageNodeAssociationSmooth, i-1)[1:5],
                                          CircShift(NodeSize, i-1)[1:5],
                                          1)
      
      if(i+4 <= length(NodeSize)){
        StageNodeAssociationSmooth[i:(i+4)] <- SmoothedAssociation
      } else {
        Tail <- i:length(StageNodeAssociationSmooth)
        Head <- 1:(5-length(Tail))
        StageNodeAssociationSmooth[Tail] <- SmoothedAssociation[1:length(Tail)]
        StageNodeAssociationSmooth[Head] <- SmoothedAssociation[(length(Tail)+1):5]
      }
      
    }
    
    
    
    ShiftMat <- NULL
    
    for(i in 1:length(StageNodeAssociationSmooth)){
      SfiftAssociation <- CircShift(StageNodeAssociationSmooth, i)
      ShiftMat <- rbind(ShiftMat, c(0, kruskal.test(1:length(StageNodeAssociationSmooth), SfiftAssociation)$p.val,
                                    aggregate(1:length(StageNodeAssociationSmooth), by=list(SfiftAssociation), mean)[,2]))
    }
    
    for(i in 1:length(StageNodeAssociationSmooth)){
      SfiftAssociation <- CircShift(rev(StageNodeAssociationSmooth), i)
      ShiftMat <- rbind(ShiftMat, c(1, kruskal.test(1:length(StageNodeAssociationSmooth), SfiftAssociation)$p.val,
                                       aggregate(1:length(StageNodeAssociationSmooth), by=list(SfiftAssociation), mean)[,2]))
    }
    
    ShiftMat <- cbind(1:length(StageNodeAssociationSmooth), ShiftMat)
    
    Candisdates <- !apply(ShiftMat[,-c(1:3)], 1, is.unsorted)
    
    if(sum(Candisdates) > 0){
      
      print("Stage V.III - Mean ordered association found")
      
      ShiftMat <- ShiftMat[Candisdates,]
      if(sum(Candisdates)>1){
        MeansSpan <- ShiftMat[,ncol(ShiftMat)] - ShiftMat[,3]
        ShiftMatSel <- ShiftMat[which.max(MeansSpan),]
      } else {
        ShiftMatSel <- ShiftMat
      }
      
      if(ShiftMatSel[2] == 1){
        UsedPath <- rev(UsedPath)
        StageNodeAssociationSmooth <- rev(StageNodeAssociationSmooth)
      }
      
      UsedPath <- CircShift(UsedPath, ShiftMatSel[1])
      StagesOnPath <- CircShift(StageNodeAssociationSmooth, ShiftMatSel[1])
      
    } else {
      
      print("Stage V.III - No mean ordered association found. Using only the 1st and last stages")
      
      MeansSpan <- ShiftMat[,ncol(ShiftMat)] - ShiftMat[,3]
      ShiftMatSel <- ShiftMat[,which.max(MeansSpan)]
      
      if(ShiftMatSel[2] == 1){
        UsedPath <- rev(UsedPath)
        StageNodeAssociationSmooth <- rev(StageNodeAssociationSmooth)
      }
      
      UsedPath <- CircShift(UsedPath, ShiftMatSel[1])
      StagesOnPath <- CircShift(StageNodeAssociationSmooth, ShiftMatSel[1])
      
    }

  } else {
    
    Optimized <- FALSE
    
    if(PathOpt == "switch"){
      print("Optimizing switch-Like behaviours")
      print("... but not now")
    }
    
    if(!Optimized){
      print("No valid optimization strategy found. Reference path will be used")
    }
    
    StagesOnPath <- rep(NA, length(UsedPath))
  }
  
  print(StagesOnPath)
  
  if(AggressiveStaging){

    print("Using aggressive staging")
    
    for(St in sort(unique(StagesOnPath))){
      StagesOnPath[intersect(which(StagesOnPath > St), 1:min(which(StagesOnPath==St)))] <- St
      StagesOnPath[intersect(which(StagesOnPath < St), max(which(StagesOnPath==St)):length(StagesOnPath))] <- St
    }
    
  }
  
  print(StagesOnPath)
  
  # Projecting on Path ---------------------------------------------------------------
  
  print("Stage VI - Path Projection")
  
  print("The following path will be used")
  print(UsedPath)
  readline("Press any key")
  
  NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
  
  NumericPath <- c(NumericPath, NumericPath[1])
  
  PathProjection <- OrderOnPath(PrinGraph = Results, Path = NumericPath, PointProjections = ProjPoints)
  
  par(mfcol=c(1,1))

  CellStages <- rep(NA, length(PathProjection$PositionOnPath))
  
  if(is.list(StageAssociation)){
    for(i in 1:length(TaxonList)){
      
      if(any(is.na(TaxonList[[i]]))){
        next()
      }
      
      CellStages[TaxonList[[i]]] <- StageAssociation$Stages[StagesOnPath[which(NumericPath == i)[1]]] 
    }
  }
  
  DF.Plot <- cbind(PathProjection$PositionOnPath/sum(PathProjection$PathLen),
                   PathProjection$DistanceFromPath,
                   CellStages)
  colnames(DF.Plot) <- c("PseudoTime", "Distance", "Stage")
  
  DF.Plot <- as.data.frame(DF.Plot)
  DF.Plot$PseudoTime <- as.numeric(as.character(DF.Plot$PseudoTime))
  DF.Plot$Distance <- as.numeric(as.character(DF.Plot$Distance))
  
  p <- ggplot2::ggplot(data = DF.Plot, mapping = aes(x = PseudoTime, y = Distance, color = Stage, label=rownames(RotatedExpression)))
  p <- p + ggplot2::geom_point() + ggplot2::geom_text(check_overlap = FALSE, hjust = "inward")
  print(p)  
  
  NodeOnGenesOnPath <- NodeOnGenes[NumericPath,]
  
  readline("Press any key")
  
  # Plot Genes ---------------------------------------------------------------
  
  print("Stage VII - Plotting")
   
  if(is.list(StageAssociation)){
    
    for (Stage in 1:length(StageAssociation$Stages)) {
      
      if(exists(paste("S", Stage, "_U", sep = ""), where=StageAssociation)) {
        
        SmoothedGenes <- GeneExpressiononPath(ExpressionData = NormExpressionMatrix, TransfData  = RotatedExpression,
                                              CellClass = CellStages, PrinGraph = Results,
                                              Projections = ProjPoints, InvTransNodes = NodeOnGenes,
                                              Genes = StageAssociation[[paste("S", Stage, "_U", sep = "")]],
                                              Path = UsedPath, Net = Net,
                                              PathType = 'Circular', Circular = TRUE,
                                              Plot = TRUE, CircExt = .1, Return.Smoother = 'loess')
        
      }
      
      if(exists(paste("S", Stage, "_D", sep = ""), where=StageAssociation)) {
        
        SmoothedGenes <- GeneExpressiononPath(ExpressionData = NormExpressionMatrix, TransfData  = RotatedExpression,
                                              CellClass = CellStages, PrinGraph = Results,
                                              Projections = ProjPoints, InvTransNodes = NodeOnGenes,
                                              Genes = StageAssociation[[paste("S", Stage, "_D", sep = "")]],
                                              Path = UsedPath, Net = Net,
                                              PathType = 'Circular', Circular = TRUE,
                                              Plot = TRUE, CircExt = .1, Return.Smoother = 'loess')
        
      }

      }
    }

  }

  
  







