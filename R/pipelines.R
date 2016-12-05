
#' Title
#'
#' @param CCList 
#' @param StageAssociation 
#' @param TopQ 
#' @param LowQ 
#' @param NodePower 
#' @param PercNorm 
#' @param MinWit 
#' @param StagingMode 
#'
#' @return
#' @export
#'
#' @examples
ReStage <- function(CCList, StageAssociation, TopQ = .9, LowQ = .1, NodePower = 0, PercNorm = TRUE,
                    MinWit = 0, StagingMode = 4, PathOpt = "Genes.PV") {
  
  # Path Selection ---------------------------------------------------------------
  
  print("Stage V - Path Selection")
  
  # Start by selecting the first available path.
  
  Pattern <- igraph::graph.ring(n = nrow(CCList$PrinGraph[[1]]$Nodes),
                                directed = FALSE, mutual = FALSE, circular = FALSE)
  
  PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = CCList$Net, graph2 = Pattern)
  
  StagingAttempts <- list()
  
  UsedPath <- PossiblePaths[[sample(1:length(PossiblePaths), 1)]]$name
  
  print(paste("Staging using the following reference path"))
  print(UsedPath)
  
  NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
  
  NodeOnGenesOnPath <- CCList$InvTransNodes[NumericPath,]
  
  Staged <- FALSE
  
  if(PathOpt == "Genes.PV" & is.list(StageAssociation)){
    
    NodeSize <- unlist(lapply(lapply(CCList$TaxonList, is.finite), sum))
    
    CutOffVar <- NULL
    
    if(exists("QVarCutOff", where=StageAssociation)){
      CutOffVar <- quantile(apply(CCList$UnscExpressionData, 2, var), as.numeric(StageAssociation$QVarCutOff))
    }
    
    StagingResults <- StagingByGenes(StageAssociation = StageAssociation, ExpressionMatrix = CCList$UnscExpressionData,
                                     NormExpressionMatrix = CCList$ExpressionData, NodeOnGenesOnPath = NodeOnGenesOnPath,
                                     UsedPath = UsedPath, NodeSize = NodeSize, NodePower = NodePower,
                                     LowQ = LowQ, TopQ = TopQ, nPoints = nrow(CCList$PrinGraph[[1]]$Nodes), MinWit = MinWit,
                                     PercNorm = PercNorm, StagingMode = StagingMode, CutOffVar = CutOffVar)
    
    # print(StagingResults)
    
    SummaryStageMat <- StagingResults$SummaryStageMat
    
    print("Optimizing path")
    
    CompactVertexStage <- NULL
    for(i in 1:length(StageAssociation$Stages)){
      CompactVertexStage <- rbind(CompactVertexStage,
                                  colSums(StagingResults$AllStg==i))
    }
    
    
    if(length(StagingResults$AllPen)>1){
      print(paste(length(StagingResults$AllPen), "minima found"))
      print("Selecting one at random")
      
      SelIdx <- sample(1:nrow(StagingResults$AllStg), 1)
      
      if(StagingResults$AllDir[SelIdx] == "Rev"){
        UsedPath <- rev(colnames(CompactVertexStage))
        StageVect <- rev(StagingResults$AllStg[SelIdx,])
      } else {
        UsedPath <- colnames(CompactVertexStage)
        StageVect <- StagingResults$AllStg[SelIdx,]
      }
    } else {
      if(StagingResults$AllDir == "Rev"){
        UsedPath <- rev(colnames(CompactVertexStage))
        StageVect <- rev(StagingResults$AllStg)
      } else {
        UsedPath <- colnames(CompactVertexStage)
        StageVect <- StagingResults$AllStg
      }
    }
    
    for (i in 1:length(StageVect)) {
      TestShift <- CircShift(StageVect, i-1)
      if(TestShift[1] == min(StageVect) & TestShift[length(TestShift)] != min(StageVect)){
        break
      }
    }
    
    UsedPath <- CircShift(UsedPath, i-1)
    StagesOnPath <- CircShift(StageVect, i-1)
    
    VertexStageMatrix <- StagingResults$AllStg[,UsedPath]
    CompactVertexStage <- CompactVertexStage[,UsedPath]
    SummaryStageMat <- SummaryStageMat[,UsedPath]
    
    print("Stage VI - Path Projection")
    
    print("The following path will be used")
    print(UsedPath)
    
    print("The following staging has been inferred")
    print(StageAssociation$Stages[StagesOnPath])
    
    NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
    
    NumericPath <- c(NumericPath, NumericPath[1])
    StagesOnPath <- c(StagesOnPath, StagesOnPath[1])
    
    Staged <- TRUE
    
  }
  
  
  if(!Staged){
    
    SummaryStageMat <- NULL
    StagesOnPath <- rep("Unassigned", nrow(CCList$PrinGraph[[1]]$Nodes))
    
    print("Staging information incompatible with current options")
    print("The defult path will be used")
    
  }
  
  
  
  # Projecting on Path ---------------------------------------------------------------
  
  PathProjection <- OrderOnPath(PrinGraph = CCList$PrinGraph[[1]], Path = NumericPath, PointProjections = CCList$Projections)
  
  # Move cells from the end to the beginning
  # Since floating point aritmetic is involved, I need to use a threshold for zero
  
  PathProjection$PositionOnPath[ abs(PathProjection$PositionOnPath - sum(PathProjection$PathLen)) < 1e-8 ] <- 0
  
  par(mfcol=c(1,1))
  
  CellOnNodes <- rep(NA, length(PathProjection$PositionOnPath))
  CellStages <- rep("Unassigned", length(PathProjection$PositionOnPath))
  GeneCount <- NULL
  CellStagesMat <- NULL
  
  NodeSize <- unlist(lapply(lapply(CCList$TaxonList, is.finite), sum))
  
  if(PathOpt == "Genes.PV" & is.list(StageAssociation)){
    for(i in 2:length(PathProjection$PathLen)){
      CellOnNodes[PathProjection$PositionOnPath >= cumsum(PathProjection$PathLen)[i-1] &
                    PathProjection$PositionOnPath < mean(cumsum(PathProjection$PathLen)[(i-1):i])] <- i-1
      CellOnNodes[PathProjection$PositionOnPath < cumsum(PathProjection$PathLen)[i] &
                    PathProjection$PositionOnPath >= mean(cumsum(PathProjection$PathLen)[(i-1):i])] <- i
    }
    
    
    CellOnNodes[CellOnNodes > nrow(CCList$PrinGraph[[1]]$Nodes)] <- CellOnNodes[CellOnNodes > nrow(CCList$PrinGraph[[1]]$Nodes)] - nrow(CCList$PrinGraph[[1]]$Nodes)
    
    CellStagesMat <- t(t(t(SummaryStageMat)/colSums(SummaryStageMat))[, CellOnNodes])
    CellStagesMat <- cbind(CellStagesMat, StagesOnPath[CellOnNodes], CCList$Grouping)
    colnames(CellStagesMat) <- c(StageAssociation$Stages, "Stage", "Group")
    
    CellStages <- StageAssociation$Stages[CellStagesMat[,length(StageAssociation$Stages)+1]]
    StagesOnPath <- StagesOnPath[-length(StagesOnPath)]
    
    
    Labels <- rownames(CCList$PCAData)
    
    DF.Plot <- cbind(PathProjection$PositionOnPath/sum(PathProjection$PathLen),
                     PathProjection$DistanceFromPath,
                     CellStages, as.character(CCList$Grouping), Labels)
    colnames(DF.Plot) <- c("PseudoTime", "Distance", "Stage", "Grouping", "Labels")
    
    DF.Plot <- as.data.frame(DF.Plot)
    DF.Plot$PseudoTime <- as.numeric(as.character(DF.Plot$PseudoTime))
    DF.Plot$Distance <- as.numeric(as.character(DF.Plot$Distance))
    
    # print(CellStages)
    # print(StageAssociation$Stages)
    
    for(Ref in rev(StageAssociation$Stages)){
      if(Ref %in% levels(DF.Plot$Stage)){
        DF.Plot$Stage <- relevel(DF.Plot$Stage, Ref)
      }
    }
    
    AtBottom <- which(DF.Plot$PseudoTime == 0)
    AtTop <- which(DF.Plot$PseudoTime == 1)
    
    print(table(DF.Plot$Grouping, DF.Plot$Stage))
    
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
    
  }

  RetList <- CCList
  
  RetList[["InferredStages"]] <- CellStages
  RetList[["Path"]] <- UsedPath
  RetList[["NumericPath"]] <- NumericPath
  RetList[["PathProjection"]] <- PathProjection
  RetList[["StagesOnPath"]] <- StagesOnPath
  RetList[["StageWitnesses"]] <- SummaryStageMat
  RetList[["StageWitnessesCount"]] <- GeneCount
  RetList[["StageWitWeigh"]] <- NodeSize^NodePower
  RetList[["CellsStagesAssociation"]] <- CellStagesMat

  return(RetList)
  
}



#' Title
#'
#' @param CCList 
#' @param nPoints 
#'
#' @return
#' @export
#'
#' @examples
ChangeNodeNumber <- function(CCList, nPoints) {
  
  # Curve Fit ---------------------------------------------------------------
  
  print("Stage V - Curve fit")
  
  Results <- computeElasticPrincipalGraph(Data = CCList$PCAData, NumNodes = nPoints, Method = 'CircleConfiguration')
  
  TaxonList <- getTaxonMap(Results[[1]], CCList$PCAData, UseR = TRUE)
  
  ProjPoints <- projectPoints(Results = Results[[1]], Data = CCList$PCAData, TaxonList=TaxonList,
                              UseR = TRUE,
                              method = 'PCALin', Dims = NULL, Debug = FALSE)
  
  Net <- ConstructGraph(Results = Results[[1]], DirectionMat = NULL)
  
  NodeOnGenes <- t(t(Results[[1]]$Nodes %*% t(CCList$PCAInfo$Comp)))
  
  RetList <- CCList
  
  RetList[["Net"]] <- Net
  RetList[["PrinGraph"]] <- Results
  RetList[["InvTransNodes"]] <- NodeOnGenes
  RetList[["TaxonList"]] <- TaxonList
  RetList[["Projections"]] <- ProjPoints
  
  RetList[["InferredStages"]] <- NULL
  RetList[["Path"]] <- NULL
  RetList[["NumericPath"]] <- NULL
  RetList[["PathProjection"]] <- NULL
  RetList[["StagesOnPath"]] <- NULL
  RetList[["StageWitnesses"]] <- NULL
  RetList[["StageWitnessesCount"]] <- NULL
  RetList[["StageWitWeigh"]] <- NULL
  RetList[["CellsStagesAssociation"]] <- NULL
  
  return(RetList)
  
}





#' Title
#'
#' @param ExpressionMatrix 
#' @param GeneDetectedFilter 
#' @param GeneCountFilter 
#' @param MinCellExp 
#' @param LogTranform 
#' @param QuantNorm 
#' @param GeneSet 
#' @param Centering 
#' @param Scaling 
#'
#' @return
#' @export
#'
#' @examples
Filter <- function(ExpressionMatrix, GeneDetectedFilter, GeneCountFilter, MinCellExp, LogTranform,
                   QuantNorm, GeneSet, Centering, Scaling, Grouping) {
  
  # Cell filtering ---------------------------------------------------------------
  
  print("Stage I - Cell filtering")
  
  OutExpr <- scater::isOutlier(rowSums(ExpressionMatrix>0), nmads = GeneDetectedFilter)
  
  OutCount <- scater::isOutlier(rowSums(ExpressionMatrix>0), nmads = GeneCountFilter)
  
  print(paste(sum(OutExpr & OutCount), "Cells will be removed due to both filtering"))
  print(paste(sum(OutExpr & !OutCount), "Cells will be removed due to gene count filtering"))
  print(paste(sum(OutCount & !OutExpr), "Cells will be removed due to read count filtering"))
  print(paste(sum(!(OutExpr & OutCount)), "Cells will be used for analysis")) 
  
  Grouping <- Grouping[!(OutExpr & OutCount)]
  NormExpressionMatrix <- ExpressionMatrix[!(OutExpr & OutCount),]
  
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
  
  UnScaledNormExpressionMatrix <- NormExpressionMatrix
  
  if(QuantNorm){
    Rname <- rownames(NormExpressionMatrix)
    Cname <- colnames(NormExpressionMatrix)
    NormExpressionMatrix <- preprocessCore::normalize.quantiles(NormExpressionMatrix)
    rownames(NormExpressionMatrix) <- Rname
    colnames(NormExpressionMatrix) <- Cname
  }
  
  if(!is.null(GeneSet)){
    SelGenes <- intersect(colnames(NormExpressionMatrix), GeneSet)
    
    print(paste(length(SelGenes), "genes selected for analysis"))
    
    if(length(SelGenes)==0){
      return()
    }
    
    NormExpressionMatrix <- NormExpressionMatrix[,SelGenes]
    
  }
  
  # Gene Scaling ---------------------------------------------------------------
  
  print("Stage III - Gene scaling")
  
  if(Centering){
    print("Gene expression will be centered")
  } 
  
  if(Scaling){
    print("Gene expression will be scaled")
  } 
  
  NormExpressionMatrix <- scale(x = NormExpressionMatrix, center = Centering, scale =  Scaling)
  
  print("Stage IV - PCA")
  
  if(is.null(nDim)){
    nDim <- min(dim(NormExpressionMatrix))
  }
  
  NormExpressionMatrixPCA <- SelectComputePCA(NormExpressionMatrix,
                                              Components = nDim, Method = 'base-svd',
                                              center = FALSE, scale. = FALSE)
  
  RotatedExpression <- NormExpressionMatrix %*% NormExpressionMatrixPCA$Comp
  
  return(list(ExpressionData = NormExpressionMatrix, UnscExpressionData = UnScaledNormExpressionMatrix,
              PCAData = RotatedExpression,
              PCAInfo = NormExpressionMatrixPCA, Grouping = Grouping, TaxonList = NULL,
              InferredStages = NULL, PrinGraph = NULL, Projections = NULL, InvTransNodes = NULL,
              Path = NULL, NumericPath = NULL, Net = NULL, PathProjection = NULL, StagesOnPath = NULL,
              StageWitnesses = NULL, StageWitnessesCount = NULL, StageWitWeigh = NULL,
              CellsStagesAssociation = NULL))
  
}




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
StudyCellCycles <- function(ExpressionMatrix, Grouping, GeneSet = NULL, QuantNorm = FALSE,
                            StageAssociation = NULL, TopQ = .9, LowQ = .1, NodePower = 0, PercNorm = TRUE,
                            MinWit = 0, StagingMode = 4,
                            PathOpt = "Genes.PV", GeneOpt = 10, 
                            GeneDetectedFilter = 2.5, GeneCountFilter = 2.5,
                            MinCellExp = 2, VarFilter = 0, LogTranform = TRUE, Centering = FALSE,
                            Scaling = FALSE, nDim = NULL, DimPrj = NULL, nPoints = 20,
                            Data.Return = FALSE, GeneToPlot = 10, ThrNumb = NULL,
                            Interactive = TRUE) {
  
  print(paste("Expression matrix contains", nrow(ExpressionMatrix), "cells and", ncol(ExpressionMatrix), "genes"))
  
  if(length(Grouping) != nrow(ExpressionMatrix)){
    stop("Grouping vector incompatible with expression matrix")
  }
  
  
  
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
  print(paste(sum(OutExpr & !OutCount), "Cells will be removed due to gene count filtering"))
  print(paste(sum(OutCount & !OutExpr), "Cells will be removed due to read count filtering"))
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
  
  UnScaledNormExpressionMatrix <- NormExpressionMatrix
  
  if(QuantNorm){
    Rname <- rownames(NormExpressionMatrix)
    Cname <- colnames(NormExpressionMatrix)
    NormExpressionMatrix <- preprocessCore::normalize.quantiles(NormExpressionMatrix)
    rownames(NormExpressionMatrix) <- Rname
    colnames(NormExpressionMatrix) <- Cname
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
  
  
  NormExpressionMatrix <- scale(x = NormExpressionMatrix, center = Centering, scale =  Scaling)
  
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
  
  TaxonList <- getTaxonMap(Results[[1]], RotatedExpression, UseR = TRUE)
  
  ProjPoints <- projectPoints(Results = Results[[1]], Data = RotatedExpression, TaxonList=TaxonList,
                              UseR = TRUE,
                              method = 'PCALin', Dims = DimPrj, Debug = FALSE)
  
  Net <- ConstructGraph(Results = Results[[1]], DirectionMat = NULL)
  
  if(Interactive){
    
    par(mfcol=c(1,2))
    
    accuracyComplexityPlot(Results[[1]], AdjFactor = 1, Mode = 'LocMin')
    plotMSDEnergyPlot(Results[[1]])
    
    
    ColCells <- rainbow(length(unique(Grouping)))
    Col = ColCells[as.integer(factor(ColCells))]
    
    par(mfcol=c(1,1))
    
    plotData2D(Data = RotatedExpression, PrintGraph = Results[[1]], Col = ColCells, NodeSizeMult = 0.1,
               Main = "All genes", Plot.ly = TRUE, GroupsLab = Grouping,
               Xlab = paste("PC1 (", signif(100*NormExpressionMatrixPCA$ExpVar[1], 4), "%)", sep=''),
               Ylab = paste("PC2 (", signif(100*NormExpressionMatrixPCA$ExpVar[2], 4), "%)", sep=''))
    
    readline("Press any key")
    
    # arrows(x0 = RotatedExpression[,1], y0 = RotatedExpression[,2],
    #        x1 = ProjPoints$PointsOnEdgesCoords[,1], y1 = ProjPoints$PointsOnEdgesCoords[,2], length = 0)
    
    PieNetInfo <- plotPieNet(Data = RotatedExpression, Results = Results[[1]], Categories = Grouping, Graph = Net,
               TaxonList = TaxonList, NodeSizeMult = 6)
    
    legend(x = "center", legend = levels(Grouping), fill = PieNetInfo$ColInfo)
    
    readline("Press any key")
  
  }
  
  # Path Selection ---------------------------------------------------------------
  
  print("Stage V - Path Selection")
  
  # Start by selecting the first available path.
  
  Pattern <- igraph::graph.ring(n = nPoints, directed = FALSE, mutual = FALSE, circular = FALSE)
  
  PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = Net, graph2 = Pattern)
  
  NodeOnGenes <- t(t(Results[[1]]$Nodes %*% t(NormExpressionMatrixPCA$Comp)))
  
  StagingAttempts <- list()
  
  UsedPath <- PossiblePaths[[sample(1:length(PossiblePaths), 1)]]$name
  
  print(paste("Staging using the following reference path"))
  print(UsedPath)
  
  NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
  
  NodeOnGenesOnPath <- NodeOnGenes[NumericPath,]
  
  Staged <- FALSE
  
  if(PathOpt == "Genes.PV" & is.list(StageAssociation)){
    
    NodeSize <- unlist(lapply(lapply(TaxonList, is.finite), sum))
    
    CutOffVar <- NULL
    
    if(exists("QVarCutOff", where=StageAssociation)){
      CutOffVar <- quantile(apply(UnScaledNormExpressionMatrix, 2, var), as.numeric(StageAssociation$QVarCutOff))
    }
    
    StagingResults <- StagingByGenes(StageAssociation = StageAssociation, ExpressionMatrix = UnScaledNormExpressionMatrix,
                                     NormExpressionMatrix = NormExpressionMatrix, NodeOnGenesOnPath = NodeOnGenesOnPath,
                                     UsedPath = UsedPath, NodeSize = NodeSize, NodePower = NodePower,
                                     LowQ = LowQ, TopQ = TopQ, nPoints = nPoints, MinWit = MinWit,
                                     PercNorm = PercNorm, StagingMode = StagingMode, CutOffVar = CutOffVar)
    
    print(StagingResults)
    
    SummaryStageMat <- StagingResults$SummaryStageMat
    
    print("Optimizing path")
    
    CompactVertexStage <- NULL
    for(i in 1:length(StageAssociation$Stages)){
      CompactVertexStage <- rbind(CompactVertexStage,
                                  colSums(StagingResults$AllStg==i))
    }
    
    if(Interactive){
      barplot(CompactVertexStage, beside = FALSE, col = rainbow(length(StageAssociation$Stages)),
              las=2, ylab = 'Number of staging attempts', main = "Techinical staging uncertainty")
      
      barplot(CompactVertexStage, beside = TRUE, col = rainbow(length(StageAssociation$Stages)),
              las=2, ylab = 'Number of staging attempts', main = "Techinical staging uncertainty")
    }
    
    
    if(length(StagingResults$AllPen)>1){
      print(paste(length(StagingResults$AllPen), "minima found"))
      print("Selecting one at random")
      
      SelIdx <- sample(1:nrow(StagingResults$AllStg), 1)
      
      if(StagingResults$AllDir[SelIdx] == "Rev"){
        UsedPath <- rev(colnames(CompactVertexStage))
        StageVect <- rev(StagingResults$AllStg[SelIdx,])
      } else {
        UsedPath <- colnames(CompactVertexStage)
        StageVect <- StagingResults$AllStg[SelIdx,]
      }
    } else {
      if(StagingResults$AllDir == "Rev"){
        UsedPath <- rev(colnames(CompactVertexStage))
        StageVect <- rev(StagingResults$AllStg)
      } else {
        UsedPath <- colnames(CompactVertexStage)
        StageVect <- StagingResults$AllStg
      }
    }
    
    for(i in 1:length(StageVect)) {
      TestShift <- CircShift(as.vector(StageVect), n = i-1)
      print("Testing")
      print(TestShift)
      print(i-1)
      if(TestShift[1] == min(StageVect) & TestShift[length(TestShift)] != min(StageVect)){
        break
      }
    }
    
    UsedPath <- CircShift(as.vector(UsedPath), i-1)
    StagesOnPath <- CircShift(as.vector(StageVect), i-1)
    
    VertexStageMatrix <- StagingResults$AllStg[,UsedPath]
    CompactVertexStage <- CompactVertexStage[,UsedPath]
    SummaryStageMat <- SummaryStageMat[,UsedPath]
    
    
    if(Interactive){
      barplot(CompactVertexStage, beside = FALSE, col = rainbow(length(StageAssociation$Stages)), las=2)
      
      barplot(CompactVertexStage, beside = TRUE, col = rainbow(length(StageAssociation$Stages)), las= 2)
    }
    
    print("Stage VI - Path Projection")
    
    print("The following path will be used")
    print(UsedPath)
    
    print("The following staging has been inferred")
    print(StageAssociation$Stages[StagesOnPath])
    
    NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
    
    NumericPath <- c(NumericPath, NumericPath[1])
    StagesOnPath <- c(StagesOnPath, StagesOnPath[1])
    
    if(Interactive){
      readline("Press any key")
    }
    
    Staged <- TRUE
    
  }
  
  
  if(!Staged){
    
    SummaryStageMat <- NULL
    StagesOnPath <- rep("Unassigned", nPoints)
    
    print("Staging information incompatible with current options")
    print("The defult path will be used")
    
  }
  
  
  
  
  
  # Projecting on Path ---------------------------------------------------------------
  
  PathProjection <- OrderOnPath(PrinGraph = Results[[1]], Path = NumericPath, PointProjections = ProjPoints)
  
  # Move cells from the end to the beginning
  # Since floating point aritmetic is involved, I need to use a threshold for zero
  
  PathProjection$PositionOnPath[ abs(PathProjection$PositionOnPath - sum(PathProjection$PathLen)) < 1e-8 ] <- 0
  
  par(mfcol=c(1,1))

  CellOnNodes <- rep(NA, length(PathProjection$PositionOnPath))
  CellStages <- rep("Unassigned", length(PathProjection$PositionOnPath))
  GeneCount <- NULL
  CellStagesMat <- NULL
  
  NodeSize <- unlist(lapply(lapply(TaxonList, is.finite), sum))
  
  if(PathOpt == "Genes.PV" & is.list(StageAssociation)){
    for(i in 2:length(PathProjection$PathLen)){
      CellOnNodes[PathProjection$PositionOnPath >= cumsum(PathProjection$PathLen)[i-1] &
                   PathProjection$PositionOnPath < mean(cumsum(PathProjection$PathLen)[(i-1):i])] <- i-1
      CellOnNodes[PathProjection$PositionOnPath < cumsum(PathProjection$PathLen)[i] &
                   PathProjection$PositionOnPath >= mean(cumsum(PathProjection$PathLen)[(i-1):i])] <- i
    }
    
    
    CellOnNodes[CellOnNodes > nPoints] <- CellOnNodes[CellOnNodes > nPoints] - nPoints
    
    CellStagesMat <- t(t(t(SummaryStageMat)/colSums(SummaryStageMat))[, CellOnNodes])
    CellStagesMat <- cbind(CellStagesMat, StagesOnPath[CellOnNodes], Grouping)
    colnames(CellStagesMat) <- c(StageAssociation$Stages, "Stage", "Group")
    
    if(Interactive){
      barplot(t(CellStagesMat[order(PathProjection$PositionOnPath),1:length(StageAssociation$Stages)]),
              horiz = TRUE, las=2, col=rainbow(length(StageAssociation$Stages)),
              names.arg = StageAssociation$Stages[CellStagesMat[order(PathProjection$PositionOnPath),length(StageAssociation$Stages)+1]])
      abline(v = 100/length(StageAssociation$Stages)/100)
    }
    
    CellStages <- StageAssociation$Stages[CellStagesMat[,length(StageAssociation$Stages)+1]]
    StagesOnPath <- StagesOnPath[-length(StagesOnPath)]
    
    Labels <- rownames(RotatedExpression)
    
    DF.Plot <- cbind(PathProjection$PositionOnPath/sum(PathProjection$PathLen),
                     PathProjection$DistanceFromPath,
                     CellStages, as.character(Grouping), Labels)
    colnames(DF.Plot) <- c("PseudoTime", "Distance", "Stage", "Grouping", "Labels")
    
    DF.Plot <- as.data.frame(DF.Plot)
    DF.Plot$PseudoTime <- as.numeric(as.character(DF.Plot$PseudoTime))
    DF.Plot$Distance <- as.numeric(as.character(DF.Plot$Distance))
    
    # print(CellStages)
    # print(StageAssociation$Stages)
    
    for(Ref in rev(StageAssociation$Stages)){
      if(Ref %in% levels(DF.Plot$Stage)){
        DF.Plot$Stage <- relevel(DF.Plot$Stage, Ref)
      }
    }
    
    # AtBottom <- which(DF.Plot$PseudoTime == 0)
    # AtTop <- which(DF.Plot$PseudoTime == 1)
    
    print(table(DF.Plot$Grouping, DF.Plot$Stage))
    
    if(Interactive & length(unique(Grouping[!is.na(Grouping)])) > 1){
      gplots::heatmap.2(table(DF.Plot$Grouping, DF.Plot$Stage))
    }
    
    # if(length(AtBottom) > 0){
    #   BottomCells <- DF.Plot[AtBottom,]
    #   BottomCells$PseudoTime <- 1
    #   DF.Plot <- rbind(DF.Plot, BottomCells)
    # }
    
    # if(length(AtTop) > 0){
    #   TopCells <- DF.Plot[AtTop,]
    #   TopCells$PseudoTime <- 0
    #   DF.Plot <- rbind(DF.Plot, TopCells)
    # }
    
    if(Interactive){
      
      p <- ggplot2::ggplot(data = DF.Plot, mapping = ggplot2::aes(x = PseudoTime, y = Distance, color = Stage, label=Labels, shape=Grouping))
      p <- p + ggplot2::geom_point(size = 5) + ggplot2::geom_text(check_overlap = FALSE, hjust = "inward")
      print(p)  
      
      NodeOnGenesOnPath <- NodeOnGenes[NumericPath,]
      
      DF2.ToPlot <- cbind(rep(1:ncol(SummaryStageMat), nrow(SummaryStageMat)),
                          StagingResults$NoNormWeigth*as.vector(t(SummaryStageMat)),
                          rep(StageAssociation$Stages, each=ncol(SummaryStageMat)),
                          rep(StageAssociation$Stages[StagesOnPath], nrow(SummaryStageMat)))
      colnames(DF2.ToPlot) <- c("Node", "Percentage", "WStage", "AStage")
      
      DF2.ToPlot <- data.frame(DF2.ToPlot)
      DF2.ToPlot$Node <- as.numeric(as.character(DF2.ToPlot$Node))
      DF2.ToPlot$Percentage <- as.numeric(as.character(DF2.ToPlot$Percentage))
      
      for(Ref in rev(StageAssociation$Stages)){
        if(Ref %in% levels(DF2.ToPlot$WStage)){
          DF2.ToPlot$WStage <- relevel(DF2.ToPlot$WStage, Ref)
        }
      }
      
      for(Ref in rev(StageAssociation$Stages)){
        if(Ref %in% levels(DF2.ToPlot$AStage)){
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
    
    
  }
  
  # Finetuning the positions depending on the staging of the 1st stage 
  
  print("Finetuning cell positions on the path")
  
  ToReposition <- which(CellOnNodes == 1 & PathProjection$PositionOnPath >= sum(PathProjection$PathLen[-nPoints]))
  
  if(length(ToReposition)>0){
    
    # There are cells that are at the end and should be at the beginning
    
    AdjPar <- min(PathProjection$PositionOnPath[ToReposition] - sum(PathProjection$PathLen))
    
    # I'm going to shift down all of the projection structure by AdjPar. Note that the sign are reversed because AdjPar is negative
    
    PathProjection$PositionOnPath[ToReposition] <- PathProjection$PositionOnPath[ToReposition] - sum(PathProjection$PathLen)
    
    PathProjection$PositionOnPath <- PathProjection$PositionOnPath - AdjPar
    PathProjection$PathLen[2] <- PathProjection$PathLen[2] - AdjPar
    PathProjection$PathLen[nPoints + 1] <- PathProjection$PathLen[nPoints + 1] + AdjPar
  }
  
  
  if(Interactive){
    readline("Press any key")
  }
  
  # Plot Genes ---------------------------------------------------------------
  
  if(Interactive){
    
    
    print("Stage VII - Plotting")
    
    if(PathOpt == "Genes.PV" & is.list(StageAssociation)){
      
      print("Stage VII.I - Plotting staging genes")
      
      for (Stage in 1:length(StageAssociation$Stages)) {
        
        if(exists(paste("S", Stage, "_U", sep = ""), where=StageAssociation)) {
          
          AvailableGenes <- intersect(StageAssociation[[paste("S", Stage, "_U", sep = "")]], colnames(NodeOnGenesOnPath))
          
          if(length(AvailableGenes)>0 & !is.null(CutOffVar)){
            AvailableGenes <- AvailableGenes[apply(NormExpressionMatrix[,AvailableGenes], 2, var) > CutOffVar]
          }
          
          SmoothedGenes <- GeneExpressiononPath(ExpressionData = NormExpressionMatrix, TransfData  = RotatedExpression,
                                                CellClass = factor(CellStages, levels = StageAssociation$Stages),
                                                PrinGraph = Results[[1]],
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
                                                PrinGraph = Results[[1]],
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
                                              PrinGraph = Results[[1]],
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
                                            PrinGraph = Results[[1]],
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
    return(list(ExpressionData = NormExpressionMatrix, PCAData = RotatedExpression, PCAInfo = NormExpressionMatrixPCA,
                Grouping = Grouping, UnscExpressionData = UnScaledNormExpressionMatrix,
                InferredStages = CellStages, PrinGraph = Results, Projections = ProjPoints, InvTransNodes = NodeOnGenes,
                Path = UsedPath, NumericPath = NumericPath, Net = Net, PathProjection = PathProjection, StagesOnPath = StagesOnPath,
                StageWitnesses = SummaryStageMat, StageWitnessesCount = GeneCount, StageWitWeigh = NodeSize^NodePower,
                CellsStagesAssociation = CellStagesMat, TaxonList = TaxonList))
    }

}




#' Title
#'
#' @param CCStruct 
#'
#' @return
#' @export
#'
#' @examples
MakeCCSummaryMatrix <- function(CCStruct) {
  
  pb <- txtProgressBar(min = 0, max = ncol(CCStruct$UnscExpressionData), initial = 0, style = 3)
  
  RunTests <- function(idx) {
    
    setTxtProgressBar(pb, idx)
    
    RowVect <- rep(NA, 15)
    
    StageCount <- table(factor(CCStruct$InferredStages, levels = c("G0", "G1", "S", "G2", "M")))
    
    if(all(StageCount > 3)){
      RowVect[1] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages))$p.val
    }
    
    if(StageCount["G0"] > 3){
      RowVect[2] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages == "G0"))$p.val
    }
    
    if(StageCount["G1"] > 3){
      RowVect[3] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages == "G1"))$p.val
    }
   
    if(StageCount["S"] > 3){
      RowVect[4] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages == "S"))$p.val
    }
    
    if(StageCount["G2"] > 3){
      RowVect[5] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages == "G2"))$p.val
    }
    
    if(StageCount["M"] > 3){
      RowVect[6] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages == "M"))$p.val
    }
    
    if(StageCount["G0"] + StageCount["G1"] > 3){
      RowVect[7] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages %in% c("G0", "G1")))$p.val
    }
    
    
    if(StageCount["G1"] + StageCount["S"] > 3){
      RowVect[8] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages %in% c("G1", "S")))$p.val
    }
    
    
    if(StageCount["S"] + StageCount["G2"] > 3){
      RowVect[9] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages %in% c("S", "G2")))$p.val
    }
    
    if(StageCount["G2"] + StageCount["M"] > 3){
      RowVect[10] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages %in% c("G2", "M")))$p.val
    }
    
    if(StageCount["M"] + StageCount["G0"] > 3){
      RowVect[11] <- kruskal.test(CCStruct$UnscExpressionData[, idx], factor(CCStruct$InferredStages %in% c("M", "G0")))$p.val
    }
    
    RowVect[12] <- var(CCStruct$UnscExpressionData[, idx])
    RowVect[13] <- mean(CCStruct$UnscExpressionData[, idx])
    
    RowVect[14] <- IQR(CCStruct$UnscExpressionData[, idx])
    RowVect[15] <- median(CCStruct$UnscExpressionData[, idx])
    
    return(RowVect)
    
  }
  
  SummaryTable <- sapply(1:ncol(CCStruct$UnscExpressionData), RunTests)
  
  SummaryTable <- t(SummaryTable)
  
  SummaryTable <- cbind(colnames(CCStruct$UnscExpressionData), SummaryTable)
  colnames(SummaryTable) <- c("Gene", "KT all", "KT G0", "KT G1", "KT S", "KT G2", "KT M", "KT G0+G1", "KT G1+S", "KT S+G2", "KT G2+M", "KT M+G0",
                              "Var", "Mean", "IQR", "Median")
  
  return(SummaryTable)
  
}



  
