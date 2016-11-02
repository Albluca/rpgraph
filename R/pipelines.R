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
                            StageAssociation = NULL, TopQ = .9, LowQ = .1, NodePower = 2, PercNorm = TRUE,
                            MinWit = 0, ReStage = 100,
                            PathOpt = "switch", GeneOpt = 10,
                            GeneDetectedFilter = 2.5, GeneCountFilter = 2.5,
                            MinCellExp = 1, VarFilter = 0, LogTranform = TRUE, Centering = FALSE,
                            Scaling = FALSE, nDim = NULL, nPoints = 20,
                            Data.Return = FALSE, GeneToPlot = 10, ThrNumb = NULL,
                            Interactive = TRUE) {
  
  print(paste("Expression matrix contains", nrow(ExpressionMatrix), "cells and", ncol(ExpressionMatrix), "genes"))
  
  # Cell filtering ---------------------------------------------------------------
  
  print("Stage I - Cell filtering")
  
  if(Interactive){
    
    par(mfcol=c(1,2))
    
    hist(apply(ExpressionMatrix > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
         freq = TRUE, ylab = "Number of cells")
  }
  
  OutExpr <- scater::isOutlier(rowSums(ExpressionMatrix>0), nmads = GeneDetectedFilter)
  
  
  if(Interactive){
    hist(apply(ExpressionMatrix, 1, sum), main = "Reads per cell", xlab = "Reads count",
         freq = TRUE, ylab = "Number of cells")
    
    par(mfcol=c(1,1))  
  }
  
  OutCount <- scater::isOutlier(rowSums(ExpressionMatrix>0), nmads = GeneCountFilter)

  
  print(paste(sum(OutExpr & OutCount), "Cells will be removed due to both filtering"))
  print(paste(sum(OutExpr), "Cells will be removed due to gene count filtering"))
  print(paste(sum(OutCount), "Cells will be removed due to read count filtering"))
  print(paste(sum(!(OutExpr & OutCount)), "Cells will be used for analysis")) 

  Grouping <- Grouping[!(OutExpr & OutCount)]
  NormExpressionMatrix <- ExpressionMatrix[!(OutExpr & OutCount),]
  
  if(Interactive){
    readline("Press any key")
    
    par(mfcol=c(1,2))
    
    hist(apply(NormExpressionMatrix > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
         freq = TRUE, ylab = "Number of cells (After filtering)")
    
    hist(apply(NormExpressionMatrix, 1, sum), main = "Reads per cell", xlab = "Reads count",
         freq = TRUE, ylab = "Number of cells (After filtering)")
  
    par(mfcol=c(1,1))
  }
  
  
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

  VarVect <- apply(NormExpressionMatrix, 2, var)
  
  if(Interactive){
    
    par(mfcol=c(1,2))
    
    plot(density(VarVect, from=min(VarVect)), xlab = "Variance", main = "Variance distribution")
    abline(v=quantile(VarVect, c(.01, .05, .25, .50)), lty=2, col=c('red', 'green', 'blue', "black"))
    legend(x = "topright", legend = c("1% Q", "5% Q", "25% Q", "50% Q"),
           col=c('red', 'green', 'blue', "black"), lty=2)
    
    plot(density(VarVect, from=min(VarVect)), xlab = "Variance (1% - 75% Q)",
         xlim = quantile(VarVect, c(.01, .75)), main = "Variance distribution")
    abline(v=quantile(VarVect, c(.01, .05, .25, .50)), lty=2, col=c('red', 'green', 'blue', "black"))
    
    legend(x = "topright", legend = c("1% Q", "5% Q", "25% Q", "50% Q"),
           col=c('red', 'green', 'blue', "black"), lty=2)
    
    par(mfcol=c(1,1))
    
  }
  
  
  if(!is.null(GeneSet)){
    SelGenes <- intersect(colnames(NormExpressionMatrix), GeneSet)
    
    print(paste(length(SelGenes), "genes selected for analysis"))
    
    if(length(SelGenes)==0){
      return()
    }
    
    NormExpressionMatrix <- NormExpressionMatrix[,SelGenes]
    
  }
  
  if(Interactive){
    readline("Press any key")
  }
  
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
                                              center = FALSE, scale. = FALSE)
  
  RotatedExpression <- NormExpressionMatrix %*% NormExpressionMatrixPCA$Comp
  
  # Curve Fit ---------------------------------------------------------------
  
  print("Stage V - Curve fit")
  
  Results <- computeElasticPrincipalGraph(Data = RotatedExpression, NumNodes = nPoints, Method = 'CircleConfiguration')
  
  TaxonList <- getTaxonMap(Results, RotatedExpression, UseR = TRUE)
  
  ProjPoints <- projectPoints(Results = Results, Data = RotatedExpression, TaxonList=TaxonList,
                              UseR = TRUE,
                              method = 'PCALin', Dims = NULL)
  
  Net <- ConstructGraph(Results = Results, DirectionMat = NULL)
  
  if(Interactive){
    
    par(mfcol=c(1,2))
    
    accuracyComplexityPlot(Results, AdjFactor = 1, Mode = 'LocMin')
    plotMSDEnergyPlot(Results)
    
    
    ColCells <- rainbow(length(unique(Grouping)))
    Col = ColCells[as.integer(factor(ColCells))]
    
    par(mfcol=c(1,1))
    
    plotData2D(Data = RotatedExpression, PrintGraph = Results, Col = ColCells, NodeSizeMult = 0.1,
               Main = "All genes", Plot.ly = TRUE, GroupsLab = Grouping,
               Xlab = paste("PC1 (", signif(100*NormExpressionMatrixPCA$ExpVar[1], 4), "%)", sep=''),
               Ylab = paste("PC2 (", signif(100*NormExpressionMatrixPCA$ExpVar[2], 4), "%)", sep=''))
    
    # arrows(x0 = RotatedExpression[,1], y0 = RotatedExpression[,2],
    #        x1 = ProjPoints$PointsOnEdgesCoords[,1], y1 = ProjPoints$PointsOnEdgesCoords[,2], length = 0)
    
    readline("Press any key")
  
  }
  
  # Path Selection ---------------------------------------------------------------
  
  print("Stage V - Path Selection")
  
  # Start be selecting the first available path.
  
  Pattern <- igraph::graph.ring(n = nPoints, directed = FALSE, mutual = FALSE, circular = FALSE)
  
  PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = Net, graph2 = Pattern)
  
  NodeOnGenes <- t(t(Results$Nodes %*% t(NormExpressionMatrixPCA$Comp)))
  
  StagingAttempts <- list()
  
  for(Count in 1:ReStage){
    
    UsedPath <- PossiblePaths[[sample(1:length(PossiblePaths), 1)]]$name
    
    print(paste("Staging", Count, "using the following reference path"))
    print(UsedPath)
    
    NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
    
    NodeOnGenesOnPath <- NodeOnGenes[NumericPath,]
    
    
    if(is.list(StageAssociation)){
      
      print("Gene/Stage information found. Trying to Optimize")
      
      StageMatU <- NULL
      StageMatD <- NULL
      GeneCount <- rep(0, length(StageAssociation$Stages))
      
      CutOffVar <- NULL
      
      if(exists("QVarCutOff", where=StageAssociation)){
        CutOffVar <- quantile(apply(NormExpressionMatrix, 2, var), as.numeric(StageAssociation$QVarCutOff))
      }
      
      print("Stage V.I - Associating peaks and valleys")
      
      for (Stage in 1:length(StageAssociation$Stages)) {
        
        if(exists(paste("S", Stage, "_U", sep = ""), where=StageAssociation)) {
          
          StageGenes <- unlist(StageAssociation[paste("S", Stage, "_U", sep = "")], use.names = FALSE)
          
          AvailableGenes <- intersect(StageGenes, colnames(NodeOnGenesOnPath))
          
          if(length(AvailableGenes)>0 & !is.null(CutOffVar)){
            RestrictedExpressionMatrix <- NormExpressionMatrix[,AvailableGenes]
            dim(RestrictedExpressionMatrix) <- c(length(RestrictedExpressionMatrix)/length(AvailableGenes),
                                                 length(AvailableGenes))
            AvailableGenes <- AvailableGenes[apply(RestrictedExpressionMatrix, 2, var) > CutOffVar]
            print(paste("S", Stage, "_U: ", length(AvailableGenes), " passed cutoff selection", sep = ""))
          }
          
          if(length(AvailableGenes)>0){
            
            StageTracks <- NodeOnGenesOnPath[, AvailableGenes]
            dim(StageTracks) <- c(length(StageTracks)/length(AvailableGenes), length(AvailableGenes))
            
            SignStageMat <- sign(t(StageTracks) - apply(StageTracks, 2, quantile, TopQ))
            
            SignStageMat[SignStageMat <= 0] <- NA
            SignStageMat[SignStageMat > 0] <- Stage
            
            StageMatU <- rbind(StageMatU, SignStageMat)
            GeneCount[Stage] <- GeneCount[Stage] + ncol(StageTracks)
          }
          
        }
        
        if(exists(paste("S", Stage, "_D", sep = ""), where=StageAssociation)) {
          
          StageGenes <- unlist(StageAssociation[paste("S", Stage, "_D", sep = "")], use.names = FALSE)
          
          AvailableGenes <- intersect(StageGenes, colnames(NodeOnGenesOnPath))
          
          if(length(AvailableGenes)>0 & !is.null(CutOffVar)){
            RestrictedExpressionMatrix <- NormExpressionMatrix[,AvailableGenes]
            dim(RestrictedExpressionMatrix) <- c(length(RestrictedExpressionMatrix)/length(AvailableGenes),
                                                 length(AvailableGenes))
            AvailableGenes <- AvailableGenes[apply(RestrictedExpressionMatrix, 2, var) > CutOffVar]
            print(paste("S", Stage, "_D: ", length(AvailableGenes), " passed cutoff selection", sep = ""))
          }
          
          if(length(AvailableGenes)>0){
            StageTracks <- NodeOnGenesOnPath[, AvailableGenes]
            dim(StageTracks) <- c(length(StageTracks)/length(AvailableGenes), length(AvailableGenes))
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
      
      rownames(SummaryStageMat) <- StageAssociation$Stages
      
      print("Stage V.II - Maximising stage association")
      
      NodeSize <- unlist(lapply(lapply(TaxonList, is.finite), sum))
      
      NoNormSummaryStageMat <- SummaryStageMat
      NoNormWeigth <- NodeSize^NodePower
      
      SummaryStageMat[SummaryStageMat<MinWit] <- 0
      
      if(PercNorm){
        SummaryStageMat <- SummaryStageMat/GeneCount
      }
      
      SummaryStageMat[is.nan(SummaryStageMat)] <- 0
      
      # print(dim(SummaryStageMat))
      
      tictoc::tic()
      print("Direct staging")
      Staging <- FitStagesCirc(StageMatrix = SummaryStageMat,
                               NodePenalty = NodeSize^NodePower)
      tictoc::toc()
      
      tictoc::tic()
      print("Reverse staging")
      StagingRev <- FitStagesCirc(StageMatrix = SummaryStageMat[, rev(1:ncol(SummaryStageMat))],
                                  NodePenalty = rev(NodeSize)^NodePower)  
      tictoc::toc()
      
      
      if(mean(Staging$Penality) != 0 & mean(StagingRev$Penality) != 0){
        PrThr <- mean(StagingRev$Penality)/(mean(Staging$Penality)+mean(StagingRev$Penality))
      } else {
        PrThr <- 0.5
      }
      
      if(runif(1) < PrThr){
        SelId <- sample(1:length(Staging$Penality), 1, FALSE, rank(1/Staging$Penality))
        SelPenality <- Staging$Penality[[SelId]]
        StagesOnNodes <- Staging$Order[[SelId]]
      } else {
        print("Path reversal")
        SelId <- sample(1:length(StagingRev$Penality), 1, FALSE, rank(1/StagingRev$Penality))
        SelPenality <- StagingRev$Penality[[SelId]]
        StagesOnNodes <- StagingRev$Order[[SelId]]
        UsedPath <- rev(UsedPath)
      }
      
      for (i in 1:length(StagesOnNodes)) {
        TestShift <- CircShift(StagesOnNodes, i-1)
        if(TestShift[1] == min(StagesOnNodes) & TestShift[length(TestShift)] != min(StagesOnNodes)){
          break
        }
      }
      
      UsedPath <- CircShift(UsedPath, i-1)
      StagesOnPath <- CircShift(StagesOnNodes, i-1)
      NoNormSummaryStageMat <- NoNormSummaryStageMat[, CircShift(1:ncol(NoNormSummaryStageMat), i-1)]
      
      StagingAttempts[[Count]] <- list(Penality = SelPenality, UsedPath = UsedPath, StagesOnPath = StagesOnPath, NoNormSummaryStageMat = NoNormSummaryStageMat)
      
    } else {
      print("No staging information available. The defult path will be used")
      print("Future updates will include behavioural reorganization")
      
      break()
    }
    
  }
  
  print("Optimizing path")
  
  VertexStageMatrix <- matrix(rep(0, length(StagingAttempts)*nPoints), ncol = nPoints)
  colnames(VertexStageMatrix) <- StagingAttempts[[1]]$UsedPath
  
  for(i in 1:length(StagingAttempts)){
    VertexStageMatrix[i,StagingAttempts[[i]]$UsedPath] <- StagingAttempts[[i]]$StagesOnPath
  }
  
  BestBet <- NULL
  
  for(i in 1:length(StageAssociation$Stages)){
    BestBet <- c(BestBet, which.max(colSums(VertexStageMatrix==i)))
  }
  
  PathScore <- rep(0, length(PossiblePaths))
  for(k in 1:length(PossiblePaths)){
    for(i in 1:(length(BestBet)-1)){
      for(j in (i+1):length(BestBet)){
        if(which(PossiblePaths[[k]]$name ==  names(BestBet)[i]) < which(PossiblePaths[[k]]$name ==  names(BestBet)[j])){
          PathScore[k] <- PathScore[k] + 1
        }
      }
    }
  }
  
  UsedPathMatrix <- NULL
  
  for(i in 1:length(StagingAttempts)){
    UsedPathMatrix <- rbind(UsedPathMatrix, StagingAttempts[[i]]$UsedPath)
  }
  
  BestPathsID <- which(PathScore == max(PathScore))
  BestPathsCounts <- rep(0, length(BestPathsID))
  
  for(i in 1:length(BestPathsID)){
    for(j in 1:nrow(UsedPathMatrix)){
      if(all(UsedPathMatrix[,j] == PossiblePaths[[BestPathsID[i]]]$name)){
        BestPathsCounts[i] <- BestPathsCounts[i] + 1
      }
    }
  }
  
 
    
  
  # Projecting on Path ---------------------------------------------------------------
  
  print("Stage VI - Path Projection")
  
  print("The following path will be used")
  print(UsedPath)
  
  print("The following staging has been inferred")
  print(StagesOnPath)
  
  if(Interactive){
    readline("Press any key")
  }
  
  NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
  
  NumericPath <- c(NumericPath, NumericPath[1])
  StagesOnPath <- c(StagesOnPath, StagesOnPath[1])
  
  PathProjection <- OrderOnPath(PrinGraph = Results, Path = NumericPath, PointProjections = ProjPoints)
  
  # Move cells from the end to the beginning
  # Since floating point aritmetic is involved, I need to use a threshold for zero
  
  PathProjection$PositionOnPath[ abs(PathProjection$PositionOnPath - sum(PathProjection$PathLen)) < 1e-9 ] <- 0
  
  par(mfcol=c(1,1))

  CellStages <- rep(NA, length(PathProjection$PositionOnPath))
  
  
  if(is.list(StageAssociation)){
    
    for(i in 2:length(PathProjection$PathLen)){
      CellStages[PathProjection$PositionOnPath >= cumsum(PathProjection$PathLen)[i-1] &
                   PathProjection$PositionOnPath < mean(cumsum(PathProjection$PathLen)[(i-1):i])] <- StageAssociation$Stages[StagesOnPath[i-1]]
      CellStages[PathProjection$PositionOnPath < cumsum(PathProjection$PathLen)[i] &
                   PathProjection$PositionOnPath >= mean(cumsum(PathProjection$PathLen)[(i-1):i])] <- StageAssociation$Stages[StagesOnPath[i]]

    }
    
  }
  
  StagesOnPath <- StagesOnPath[-length(StagesOnPath)]
  
  Labels <- rownames(RotatedExpression)
  
  DF.Plot <- cbind(PathProjection$PositionOnPath/sum(PathProjection$PathLen),
                   PathProjection$DistanceFromPath,
                   CellStages, Grouping, Labels)
  colnames(DF.Plot) <- c("PseudoTime", "Distance", "Stage", "Grouping", "Labels")
  
  DF.Plot <- as.data.frame(DF.Plot)
  DF.Plot$PseudoTime <- as.numeric(as.character(DF.Plot$PseudoTime))
  DF.Plot$Distance <- as.numeric(as.character(DF.Plot$Distance))
  
  # print(CellStages)
  # print(StageAssociation$Stages)
  
  if(is.list(StageAssociation)){
    for(Ref in rev(StageAssociation$Stages)){
      if(Ref %in% levels(DF.Plot$Stage)){
        DF.Plot$Stage <- relevel(DF.Plot$Stage, Ref)
      }
    }
  }

  AtBottom <- which(DF.Plot$PseudoTime == 0)
  AtTop <- which(DF.Plot$PseudoTime == 1)
  
  print(table(DF.Plot$Grouping, DF.Plot$Stage))
  
  if(Interactive & length(unique(Grouping[!is.na(Grouping)])) > 1){
    gplots::heatmap.2(table(DF.Plot$Grouping, DF.Plot$Stage))
  }
  
  if(length(AtBottom) > 0){
    BottomCells <- DF.Plot[AtBottom,]
    BottomCells$PseudoTime <- 1
    DF.Plot <- rbind(DF.Plot, BottomCells)
  }
  
  if(length(AtTop) > 0){
    TopCells <- DF.Plot[AtTop,]
    TopCells$PseudoTime <- 0
    DF.Plot <- rbind(DF.Plot, TopCells)
  }

  
  if(Interactive){
    
    p <- ggplot2::ggplot(data = DF.Plot, mapping = ggplot2::aes(x = PseudoTime, y = Distance, color = Stage, label=Labels, shape=Grouping))
    p <- p + ggplot2::geom_point(size = 5) + ggplot2::geom_text(check_overlap = FALSE, hjust = "inward")
    print(p)  
    
    NodeOnGenesOnPath <- NodeOnGenes[NumericPath,]
    
    DF2.ToPlot <- cbind(rep(1:ncol(SummaryStageMat), nrow(SummaryStageMat)),
                        NoNormWeigth*as.vector(t(SummaryStageMat)),
                       rep(StageAssociation$Stages, each=ncol(SummaryStageMat)),
                       rep(StageAssociation$Stages[StagesOnPath], nrow(SummaryStageMat)))
    colnames(DF2.ToPlot) <- c("Node", "Percentage", "WStage", "AStage")
    
    DF2.ToPlot <- data.frame(DF2.ToPlot)
    DF2.ToPlot$Node <- as.numeric(as.character(DF2.ToPlot$Node))
    DF2.ToPlot$Percentage <- as.numeric(as.character(DF2.ToPlot$Percentage))
    
    for(Ref in rev(StageAssociation$Stages)){
      if(Ref %in% levels(DF2.ToPlot$WStage)){
        DF2.ToPlot$WStage <- relevel(DF2.ToPlot$WStage, Ref)
        DF2.ToPlot$AStage <- relevel(DF2.ToPlot$AStage, Ref)
      }
    }
    
    p <- ggplot2::ggplot(data = DF2.ToPlot, mapping = ggplot2::aes(x = Node, y = Percentage, shape=AStage))
    p <- p + ggplot2::geom_point(size = 2) + ggplot2::geom_line() + ggplot2::facet_wrap(~WStage)
    print(p)
    
    
    # barplot(NoNormSummaryStageMat/GeneCount, beside = TRUE)
    
    par(mfcol=c(1,2))
    
    barplot(table(DF.Plot$Stage, DF.Plot$Grouping), beside = TRUE, legend.text = levels(DF.Plot$Stage),
            main = "Stage / Phase association", ylab = "number of cells")
    
    barplot(unlist(lapply(split(PathProjection$PathLen[-1], StagesOnPath), sum)), names.arg = StageAssociation$Stages,
            ylab = "Pseudo Duration", main = "Genetic diversity of cell stages")
    
    par(mfcol=c(1,1))
    
  }
  
  
  
  
  
  if(Interactive){
    readline("Press any key")
  }
  
  # Plot Genes ---------------------------------------------------------------
  
  if(Interactive){
    
    
    print("Stage VII - Plotting")
    
    if(is.list(StageAssociation)){
      
      print("Stage VII.I - Plotting staging genes")
      
      for (Stage in 1:length(StageAssociation$Stages)) {
        
        if(exists(paste("S", Stage, "_U", sep = ""), where=StageAssociation)) {
          
          AvailableGenes <- intersect(StageAssociation[[paste("S", Stage, "_U", sep = "")]], colnames(NodeOnGenesOnPath))
          
          if(length(AvailableGenes)>0 & !is.null(CutOffVar)){
            AvailableGenes <- AvailableGenes[apply(NormExpressionMatrix[,AvailableGenes], 2, var) > CutOffVar]
          }
          
          SmoothedGenes <- GeneExpressiononPath(ExpressionData = NormExpressionMatrix, TransfData  = RotatedExpression,
                                                CellClass = factor(CellStages, levels = StageAssociation$Stages),
                                                PrinGraph = Results,
                                                Projections = ProjPoints, InvTransNodes = NodeOnGenes,
                                                Genes = AvailableGenes,
                                                Path = UsedPath, Net = Net, Title = paste(StageAssociation$Stages[Stage], "UP"),
                                                PathType = 'Circular', Circular = TRUE,
                                                Plot = TRUE, CircExt = .1, Return.Smoother = '')
          
        }
        
        if(exists(paste("S", Stage, "_D", sep = ""), where=StageAssociation)) {
          
          AvailableGenes <- intersect(StageAssociation[[paste("S", Stage, "_D", sep = "")]], colnames(NodeOnGenesOnPath))
          
          if(length(AvailableGenes)>0 & !is.null(CutOffVar)){
            AvailableGenes <- AvailableGenes[apply(NormExpressionMatrix[,AvailableGenes], 2, var) > CutOffVar]
          }
          
          SmoothedGenes <- GeneExpressiononPath(ExpressionData = NormExpressionMatrix, TransfData  = RotatedExpression,
                                                CellClass = factor(CellStages, levels = StageAssociation$Stages),
                                                PrinGraph = Results,
                                                Projections = ProjPoints, InvTransNodes = NodeOnGenes,
                                                Genes = AvailableGenes,
                                                Path = UsedPath, Net = Net, Title = paste(StageAssociation$Stages[Stage], "DOWN"),
                                                PathType = 'Circular', Circular = TRUE,
                                                Plot = TRUE, CircExt = .1, Return.Smoother = '')
          
        }
        
        if(Interactive){
          readline("Press any key")
        }
        
      }
    }
    
    print("Stage VII.II - Plotting target genes")
    
    SelectedGenes <- NULL
    
    if(length(GeneToPlot) == 1 & is.numeric(GeneToPlot)){
      
      GeneToPlot <- round(GeneToPlot)
      
      print(paste("Plotting the top", GeneToPlot, "genes with the largest variance"))
      
      VarVect <- apply(NormExpressionMatrix, 2, var)
      
      IdentifiedGenes <- names(VarVect)[order(VarVect, decreasing = TRUE)[1:GeneToPlot]]
      
      if(is.null(GeneSet)){
        
        SmoothedGenes <- GeneExpressiononPath(ExpressionData = NormExpressionMatrix, TransfData  = RotatedExpression,
                                              CellClass = factor(CellStages, levels = StageAssociation$Stages),
                                              PrinGraph = Results,
                                              Projections = ProjPoints, InvTransNodes = NodeOnGenes,
                                              Genes = IdentifiedGenes,
                                              Path = UsedPath, Net = Net, Title = "Most varying genes",
                                              PathType = 'Circular', Circular = TRUE,
                                              Plot = TRUE, CircExt = .1, Return.Smoother = '')
        
      } else {
        
        
        
        DF.Plot <- cbind(rep(IdentifiedGenes, each = nrow(NormExpressionMatrix)),
                         rep(PathProjection$PositionOnPath/sum(PathProjection$PathLen), length(IdentifiedGenes)),
                        as.vector(NormExpressionMatrix[,IdentifiedGenes]),
                        rep(CellStages, length(IdentifiedGenes)),
                        rep(Grouping, length(IdentifiedGenes)),
                        rep(Labels, length(IdentifiedGenes)))
        
        colnames(DF.Plot) <- c("Gene", "PseudoTime", "Expression", "Stage", "Grouping", "Labels")
        
        DF.Plot <- as.data.frame(DF.Plot)
        DF.Plot$PseudoTime <- as.numeric(as.character(DF.Plot$PseudoTime))
        DF.Plot$Expression <- as.numeric(as.character(DF.Plot$Expression))
        
        if(is.list(StageAssociation)){
          for(Ref in rev(StageAssociation$Stages)){
            if(Ref %in% levels(DF.Plot$Stage)){
              DF.Plot$Stage <- relevel(DF.Plot$Stage, Ref)
            }
          }
        }
        
        
        p <- ggplot2::ggplot(DF.Plot, ggplot2::aes(x = Stage, y = Expression, fill = Stage)) +
          ggplot2::geom_boxplot() + ggplot2::facet_wrap(~ Gene) +
          ggplot2::labs(title = "Most varying genes", x = "Pseudo time", y = "Gene expression")
        print(p)
        
        
        if(length(AtBottom) > 0){
          BottomCells <- DF.Plot[AtBottom,]
          BottomCells$PseudoTime <- 1
          DF.Plot <- rbind(DF.Plot, BottomCells)
        }
        
        if(length(AtTop) > 0){
          TopCells <- DF.Plot[AtTop,]
          TopCells$PseudoTime <- 0
          DF.Plot <- rbind(DF.Plot, TopCells)
        }
        
        
        p <- ggplot2::ggplot(DF.Plot, ggplot2::aes(x = PseudoTime, y = Expression)) +
          ggplot2::geom_smooth(color="black") +
          ggplot2::geom_point(mapping = ggplot2::aes(color = Stage, shape = Grouping)) +
          ggplot2::labs(title = "Most varying genes", x = "Pseudo time", y = "Gene expression") +
          ggplot2::facet_wrap( ~ Gene)
        print(p)
        
      }
  
      
    }
    
    if(is.character(GeneToPlot)){
      
      print(paste("Plotting selected genes"))
      
      SmoothedGenes <- GeneExpressiononPath(ExpressionData = NormExpressionMatrix, TransfData  = RotatedExpression,
                                            CellClass = factor(CellStages, levels = StageAssociation$Stages),
                                            PrinGraph = Results,
                                            Projections = ProjPoints, InvTransNodes = NodeOnGenes,
                                            Genes = GeneToPlot,
                                            Path = UsedPath, Net = Net, Title = "Selected genes",
                                            PathType = 'Circular', Circular = TRUE,
                                            Plot = TRUE, CircExt = .1, Return.Smoother = '')
      
      SelectedGenes <- IdentifiedGenes
      if(Interactive){
        readline("Press any key")
      }
    }
    
  }

   
  if(Interactive){
    
    
    # Look at correlations ---------------------------------------------------------------
    
    if(!is.null(ThrNumb)){
      
      print("Stage VIII - Exhamining Correlations of selected genes")
      
      tictoc::tic()
      CorrMat <- cor(NodeOnGenes)
      tictoc::toc()
      
      for (Targ in SelectedGenes) {
        
        BaseCor <- CorrMat[, Targ]
        # hist(BaseCor, xlab = "Correlation", main = paste("Correlations of", Targ))
        Selected <- order(BaseCor, decreasing = TRUE)[1:(ThrNumb+1)]
        Selected <- setdiff(Selected, which(names(BaseCor) == Targ))
        
        
        ToPlotMat <- NodeOnGenesOnPath[,Selected]
        ToPlotMat <- ToPlotMat[-nrow(ToPlotMat),]
        
        Df.ToPlot <- cbind(as.vector(ToPlotMat),
                           rep(cumsum(PathProjection$PathLen[-length(PathProjection$PathLen)])/sum(PathProjection$PathLen[-length(PathProjection$PathLen)]),
                               ncol(ToPlotMat)),
                           rep(StageAssociation$Stages[StagesOnPath], ncol(ToPlotMat)),
                           rep(colnames(ToPlotMat), each = nrow(ToPlotMat))
        )  
        
        colnames(Df.ToPlot) <- c("Expression", "PseudoTime", "Stage", "Gene")
        
        Df.ToPlot <- as.data.frame(Df.ToPlot)
        Df.ToPlot$PseudoTime <- as.numeric(as.character(Df.ToPlot$PseudoTime))
        Df.ToPlot$Expression <- as.numeric(as.character(Df.ToPlot$Expression))
        
        p <- ggplot2::ggplot(data = Df.ToPlot, mapping = ggplot2::aes(x = PseudoTime, y = Expression, color = Stage))
        p <- p + ggplot2::geom_point() + ggplot2::labs(title = paste("Genes behaving like", Targ)) + ggplot2::facet_wrap( ~ Gene)
        print(p)
        
        
        
        
        BaseCor <- CorrMat[, Targ]
        # hist(BaseCor, xlab = "Correlation", main = paste("Correlations of", Targ))
        Selected <- order(BaseCor, decreasing = FALSE)[1:ThrNumb]
        Selected <- setdiff(Selected, which(names(BaseCor) == Targ))
        
        
        ToPlotMat <- NodeOnGenesOnPath[,Selected]
        ToPlotMat <- ToPlotMat[-nrow(ToPlotMat),]
        
        Df.ToPlot <- cbind(as.vector(ToPlotMat),
                           rep(cumsum(PathProjection$PathLen[-length(PathProjection$PathLen)])/sum(PathProjection$PathLen[-length(PathProjection$PathLen)]),
                               ncol(ToPlotMat)),
                           rep(StageAssociation$Stages[StagesOnPath], ncol(ToPlotMat)),
                           rep(colnames(ToPlotMat), each = nrow(ToPlotMat))
        )  
        
        colnames(Df.ToPlot) <- c("Expression", "PseudoTime", "Stage", "Gene")
        
        Df.ToPlot <- as.data.frame(Df.ToPlot)
        Df.ToPlot$PseudoTime <- as.numeric(as.character(Df.ToPlot$PseudoTime))
        Df.ToPlot$Expression <- as.numeric(as.character(Df.ToPlot$Expression))
        
        p <- ggplot2::ggplot(data = Df.ToPlot, mapping = ggplot2::aes(x = PseudoTime, y = Expression, color = Stage))
        p <- p + ggplot2::geom_point() + ggplot2::labs(title = paste("Genes behaving opposite to", Targ)) + ggplot2::facet_wrap( ~ Gene)
        print(p)
        
        if(Interactive){
          readline("Press any key")
        }
        
      }
      
    }
    
  }
  
  if (Data.Return) {
    
    return(list(ExpressionData = NormExpressionMatrix, PCAData = RotatedExpression, Grouping = Grouping,
                InferredStages = CellStages, PrinGraph = Results, Projections = ProjPoints, InvTransNodes = NodeOnGenes,
                Path = UsedPath, NumericPath = NumericPath, Net = Net, PathProjection = PathProjection, StagesOnPath = StagesOnPath,
                StageWitnesses = SummaryStageMat, StageWitnessesCount = GeneCount, StageWitWeigh = NodeSize^NodePower))
    }

    
}
  
  







