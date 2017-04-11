################################################################################
#
# Accessory Cell Cycle Staging fucntions to alter previous stagings ------------------------------------------------------
# The functions a still under development and currently not exported
#
################################################################################

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
                                     LowQ = LowQ, TopQ = TopQ, MaxExp = MaxExp, MinExp = MinExp,
                                     nPoints = nrow(CCList$PrinGraph[[1]]$Nodes), MinWit = MinWit,
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
  
  OutCount <- scater::isOutlier(rowSums(ExpressionMatrix), nmads = GeneCountFilter)
  
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



















################################################################################
#
# Main Cell Cycle Staging ------------------------------------------------------
#
################################################################################

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
                            StageAssociation = NULL, TopQ = .9, LowQ = .1,
                            MaxExp = NULL, MinExp = NULL,
                            NodePower = 0, PercNorm = TRUE,
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
  
  
  if(is.null(MaxExp)){
    MaxExp <- quantile(as.vector(NormExpressionMatrix), .75)
  }
  
  if(is.null(MinExp)){
    MinExp <- quantile(as.vector(NormExpressionMatrix), .25)
  }
  
  
  if(LogTranform){
    print("Transforming using pseudo counts (Log10(x+1))")
    NormExpressionMatrix <- log10(NormExpressionMatrix + 1)
    MaxExp <- log10(MaxExp+1)
    MinExp <- log10(MinExp+1)
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
                                     LowQ = LowQ, TopQ = TopQ, MaxExp = MaxExp, MinExp = MinExp,
                                     nPoints = nPoints, MinWit = MinWit,
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















################################################################################
#
# Function used to project cells on a circle or lasso ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param DataSet 
#' @param GeneSet 
#' @param OutThr 
#' @param VarThr 
#' @param nNodes 
#' @param Log 
#' @param Categories 
#' @param Filter 
#' @param GraphType 
#' @param PlanVarLimit 
#' @param PlanVarLimitIC 
#' @param MinBranDiff 
#' @param InitStructNodes 
#' @param ForceLasso 
#' @param DipPVThr 
#' @param MinProlCells 
#'
#' @return
#' @export
#'
#' @examples
#' 
ProjectAndCompute <- function(DataSet, GeneSet = NULL, VarThr, nNodes, Log = TRUE, Categories = NULL,
                              Filter = TRUE, OutThr = 5, PCAFilter = FALSE, OutThrPCA = 5,
                              GraphType = 'Lasso', PlanVarLimit = .9,
                              PlanVarLimitIC = NULL, MinBranDiff = 2, InitStructNodes = 15,
                              ForceLasso = FALSE, EstProlif = "MeanPerc", QuaThr = .5,
                              NonG0Cell = NULL, DipPVThr = 1e-3, MinProlCells = 20, PCACenter = TRUE, PCAProjCenter = FALSE,
                              PlotDebug = FALSE, PlotIntermediate = FALSE){
  
  if(is.null(PlanVarLimitIC)){
    PlanVarLimitIC <- PlanVarLimit + .5*(1-PlanVarLimit)
  }
  
  DataMat <- t(DataSet)
  
  if(is.null(Categories)){
    Categories <- factor(rep("NoG", nrow(DataMat)))
  }
  
  if(!is.null(GeneSet)){
    DataMat <- DataMat[, colnames(DataMat) %in% GeneSet]
  }
  
  DataMat <- DataMat[, apply(DataMat > 0, 2, sum) > 3]
  
  
  if(PlotDebug){
    hist(apply(DataMat > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
         freq = TRUE, ylab = "Number of cells")
    
    hist(apply(DataMat, 1, sum), main = "Reads per cell", xlab = "Reads count",
         freq = TRUE, ylab = "Number of cells")
    
    hist(apply(DataMat>0, 2, sum), main = "Transcripts per cell", xlab = "Reads count",
         freq = TRUE, ylab = "Number of genes identified")
  }

  
  if(Filter){
    OutExpr <- scater::isOutlier(rowSums(DataMat>0), nmads = OutThr)
    OutCount <- scater::isOutlier(rowSums(DataMat), nmads = OutThr)
  } else {
    OutExpr <- rep(FALSE, nrow(DataMat))
    OutCount <- rep(FALSE, nrow(DataMat))
  }
  
  if(Log){
    DataMat <- log10(DataMat[!(OutExpr & OutCount),] + 1)
  } else {
    DataMat <- DataMat[!(OutExpr & OutCount),]
  }
  
  Categories <- Categories[!(OutExpr & OutCount)]
  
  if(PCAFilter){
    PCAData <- prcomp(DataMat, retx = TRUE, center = PCACenter, scale. = FALSE)
    Centroid <- colMeans(PCAData$x)
    
    Dists <- as.matrix(dist(rbind(Centroid, PCAData$x)))
    DistFromCent <- Dists[1,-1]
    
    PCAFil <- scater::isOutlier(DistFromCent, nmads = OutThrPCA)
  } else {
    PCAFil <- rep(FALSE, nrow(DataMat))
  }
  
  DataMat <- DataMat[!PCAFil, ]
  Categories <- Categories[!PCAFil]
  
  print("Estimating most likely proliferative cells")
  RankedData <- apply(DataMat, 1, rank)
  
  
  if(is.null(NonG0Cell)){
    
    if(EstProlif == "Quantile"){
      print("Using quantile separation")
      NonG0Cell <-  rownames(DataMat)[apply(DataMat, 1, median) > quantile(DataMat, QuaThr)]
    }
    
    if(EstProlif == "PercQuant"){
      print("Using quantile ordering")
      NonG0Cell <-  rownames(DataMat)[order(apply(DataMat, 1, quantile, QuaThr), decreasing = TRUE)[1:round(nrow(DataMat)/10)]]
    }
    
    if(EstProlif == "KmeansPerc"){
      print("Using kmeans on quantile data")
      Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)
      
      KM <- kmeans(x = Cellvect, centers = range(Cellvect))
      if(PlotDebug){
        boxplot(apply(DataMat/rowSums(DataMat), 1, median) ~ KM$cluster)
      }
      
      if(KM$centers[1] < KM$centers[2]){
        NonG0Cell <-  rownames(DataMat)[KM$cluster == 2]
      } else {
        NonG0Cell <-  rownames(DataMat)[KM$cluster == 1]
      }
      
    }
    
    if(EstProlif == "MeanPerc"){
      print("Using mean of quantiles")
      Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)
      NonG0Cell <- rownames(DataMat)[Cellvect > mean(Cellvect)]
      boxplot(Cellvect)
      abline(h=mean(Cellvect))
    }
    
    if(EstProlif == "AdaptiveMeanPerc"){
      print("Using adaptive mean of quantiles")
      Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr)
      AdjFact = 0
      while((sum(Cellvect != 0) < MinProlCells) & (QuaThr + AdjFact < 1)){
        AdjFact <- AdjFact + .01
        Cellvect <- apply(DataMat/rowSums(DataMat), 1, quantile, QuaThr + AdjFact)
      }
      boxplot(Cellvect, main=paste("Q =", QuaThr + AdjFact))
      abline(h=mean(Cellvect))
      NonG0Cell <- rownames(DataMat)[Cellvect > mean(Cellvect)]
    }
    
  } else {
    NonG0Cell <- intersect(NonG0Cell, rownames(DataMat))
  }
  
  print(paste(length(NonG0Cell), "strongly proliferative cells inferred"))
  
  G0Cell <-  setdiff(rownames(DataMat), NonG0Cell)

  
  
  
  
  
  if(!is.null(NonG0Cell)){
    TB <- rbind(table(Categories[rownames(DataMat) %in% NonG0Cell]), table(Categories[rownames(DataMat) %in% G0Cell]))
    barplot(TB/rowSums(TB), ylab = "Percentage of cells", xlab = "Category", beside = TRUE,
            legend.text = c("Strongly Proliferative", "Not strongly proliferative"))
    barplot(TB, ylab = "Number of cells", xlab = "Category", beside = TRUE,
            legend.text = c("Strongly Proliferative", "Not strongly proliferative"))
    
  }
  
  
  
  if(length(NonG0Cell)>MinProlCells){
    
    print("Using strongly proliferative cells for initial fitting")
    UseTree <- TRUE
    
  } else {
    
    print("Unable to find a sufficent number of strongly proliferative cells")
    
    NonG0Cell <- rownames(DataMat)
    UseTree <- FALSE
    
  }
  
  print("Transforming Data")
  
  PCAData <- prcomp(DataMat, retx = TRUE, center = PCACenter, scale.=FALSE)
  ExpVar <- PCAData$sdev^2/sum(PCAData$sdev^2)
  
  if(VarThr<1){
    nDims <- max(min(which(cumsum(ExpVar) > VarThr)), 2)
  } else {
    nDims <- max(dim(PCAData$x))
  }
  
  Data <- PCAData$x[, 1:nDims]
  
  
  if(GraphType == 'Lasso') {
    
    print("Lasso fitting")
    print("Fitting initial circle")
    
    # Step I - Construct the base circle
    
    BasicCircData <- computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,], NumNodes = 4,
                                                  Method = 'CircleConfiguration', NodeStep = 1)
    
    PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
    PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
    
    UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes) - 1
    
    while(UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC){
      
      print("Expanding initial circle")
      
      # Contiune to add node untill the circle remains planar
      
      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,],
                                                                          NumNodes = UsedNodes + 1,
                                                                          Method = 'CircleConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))
      
      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
      
      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
      
    }
    
    if(UseTree | ForceLasso){
      
      print("Branching initial circle")
      
      # Step II - Construct the initial tree
      
      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
                                                                          NumNodes = UsedNodes + 1,
                                                                          Method = 'DefaultPrincipalTreeConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))
      
     
      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
      
      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes) - 1
      
      while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){
        
        print("Keep Branching")
        
        # Step IIa - keep constructing trees
        
        BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
                                                                            NumNodes = UsedNodes + 1,
                                                                            Method = 'DefaultPrincipalTreeConfiguration',
                                                                            NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                            Edges = BasicCircData[[length(BasicCircData)]]$Edges))
        
        
        UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
        
        PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
        PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
        
        Net <- ConstructGraph(Results = BasicCircData[[length(BasicCircData)]], DirectionMat = NULL, Thr = NULL)
        if(max(igraph::degree(Net))>2){
          break()
        }
      }
      
    }
    
    # Step III - Using curves
    
    while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){
      
      print("Extending circle and branches")
      
      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
                                                                          NumNodes = nrow(BasicCircData[[length(BasicCircData)]]$Nodes)+1,
                                                                          Method = 'CurveConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))
      
      Net <- ConstructGraph(Results = BasicCircData[[length(BasicCircData)]], DirectionMat = NULL, Thr = NULL)
      
      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
      
      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
      
      if(sum(igraph::degree(Net)>2)>1){
        
        print("Multiple branches detected ... trying to select one")
        
        # There are two branhces, this is no good ...
        # Look for the biggest circle
        
        CircSize <- 3
        
        while(CircSize < igraph::vcount(Net)){
          
          Template <- igraph::graph.ring(n = CircSize, directed = FALSE, mutual = FALSE, circular = TRUE)
          
          CircSize <- CircSize + 1
          
          CircleInData <- igraph::graph.get.subisomorphisms.vf2(Net, Template)
          
          if(length(CircleInData)>0){
            # We found the cicle
            
            EndPoint <- igraph::V(Net)[igraph::degree(Net) == 1]
            Branches <- igraph::V(Net)[igraph::degree(Net)>2]
            DistMat <- igraph::distances(Net, EndPoint, Branches)
            
            BranchesLen <- apply(DistMat, 1, min)
            
            VertToRemove <- NULL
            
            if(max(BranchesLen[-which.max(BranchesLen)] - max(BranchesLen)) <= -MinBranDiff){
              
              print("Dominating brach found ... pruning the shortest one")
              
              # There is a "dominating"" branch
              ToKeep <- names(which.max(BranchesLen))
              for(i in 1:nrow(DistMat)){
                
                Source <- rownames(DistMat)[i]
                if(ToKeep == Source){
                  next()
                }
                
                Target <- names(which.min(DistMat[i,]))
                VertToRemove <- c(VertToRemove, names(igraph::get.shortest.paths(Net, from = Target, to = Source)$vpath[[1]][-1]))
              }
              
              PrunedStruct <- BasicCircData[[length(BasicCircData)]]
              
              NodesToRem <- as.integer(unlist(lapply(strsplit(VertToRemove, "V_"), "[[", 2)))
              PrunedStruct$Nodes <- PrunedStruct$Nodes[-NodesToRem, ]
              PrunedStruct$Edges <-
                PrunedStruct$Edges[!(PrunedStruct$Edges[,1] %in% NodesToRem | PrunedStruct$Edges[,2] %in% NodesToRem), ]
              
              for(VerVal in 1:max(PrunedStruct$Edges)){
                if(any(PrunedStruct$Edges == VerVal)){
                  next()
                } else {
                  PrunedStruct$Edges[PrunedStruct$Edges>VerVal] <-
                    PrunedStruct$Edges[PrunedStruct$Edges>VerVal] - 1
                }
                
              }
              
              PrunedStruct$Method <- "ManualPruning"
              
              BasicCircData[[length(BasicCircData) + 1]] <- PrunedStruct
              
              PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
              PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
              
              UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
              
              break()
              
            }
            
            
          }
          
          
        }
        
      }
      
    }
    
    FitData <- BasicCircData
  }
  
  
  
  
  
  
  
  if(GraphType == 'Circle') {
    
    print("Circle fitting")
    print("Fitting initial circle")
    
    # Step I - Construct the base circle
    
    BasicCircData <- computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,], NumNodes = 4,
                                                  Method = 'CircleConfiguration', NodeStep = 1)
    
    PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
    PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
    
    UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes) - 1
    
    while(UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC){
      
      print("Expanding initial circle")
      
      # Contiune to add node untill the circle remains planar
      
      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,],
                                                                          NumNodes = UsedNodes + 1,
                                                                          Method = 'CircleConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))
      
      
      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
      
      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
      
      print(paste("Initial circle: Nodes = ", UsedNodes, "PercPlan = ", PlanPerc))
      
    }
    
    # Step II - Using curves
    
    while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){
      
      print("Extending circle")
      
      BasicCircData <- append(BasicCircData, computeElasticPrincipalGraph(Data = Data,
                                                                          NumNodes = nrow(BasicCircData[[length(BasicCircData)]]$Nodes)+1,
                                                                          Method = 'CurveConfiguration',
                                                                          NodesPositions = BasicCircData[[length(BasicCircData)]]$Nodes,
                                                                          Edges = BasicCircData[[length(BasicCircData)]]$Edges))
      
      
      Net <- ConstructGraph(Results = BasicCircData[[length(BasicCircData)]], DirectionMat = NULL, Thr = NULL)
      
      PCAStruct <- prcomp(BasicCircData[[length(BasicCircData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
      
      UsedNodes <- nrow(BasicCircData[[length(BasicCircData)]]$Nodes)
      
      print(paste("Final circle: Nodes = ", UsedNodes, "PercPlan = ", PlanPerc))
      
    }
    
    FitData <- BasicCircData
  }
  
  
  
  
  
  
  
  
  if(GraphType == 'Line') {
    
    print("Line fitting")
    print("Fitting initial line")
    
    # Step I - Construct the base circle
    
    BasicLineData <- computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,], NumNodes = 4,
                                                  Method = 'CurveConfiguration', NodeStep = 1)
    
    PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
    PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
    
    UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes) - 1
    
    while(UsedNodes < InitStructNodes & PlanPerc > PlanVarLimitIC){
      
      print("Expanding initial line")
      
      # Contiune to add node untill the circle remains planar
      
      BasicLineData <- append(BasicLineData, computeElasticPrincipalGraph(Data = Data[rownames(Data) %in% NonG0Cell,],
                                                                          NumNodes = UsedNodes + 1,
                                                                          Method = 'CircleConfiguration',
                                                                          NodesPositions = BasicLineData[[length(BasicLineData)]]$Nodes,
                                                                          Edges = BasicLineData[[length(BasicLineData)]]$Edges))
      
      
      UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes)
      
      PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
      
    }
    
    # Step II - Using curves
    
    while(UsedNodes < nNodes & PlanPerc > PlanVarLimit){
      
      print("Extending line")
      
      BasicLineData <- append(BasicLineData, computeElasticPrincipalGraph(Data = Data,
                                                                          NumNodes = nrow(BasicLineData[[length(BasicLineData)]]$Nodes)+1,
                                                                          Method = 'CurveConfiguration',
                                                                          NodesPositions = BasicLineData[[length(BasicLineData)]]$Nodes,
                                                                          Edges = BasicLineData[[length(BasicLineData)]]$Edges))
      
     
      Net <- ConstructGraph(Results = BasicLineData[[length(BasicLineData)]], DirectionMat = NULL, Thr = NULL)
      
      PCAStruct <- prcomp(BasicLineData[[length(BasicLineData)]]$Nodes, center = TRUE, scale. = FALSE, retx = TRUE)
      PlanPerc <- sum(PCAStruct$sdev[1:2]^2)/sum(PCAStruct$sdev^2)
      
      UsedNodes <- nrow(BasicLineData[[length(BasicLineData)]]$Nodes)
      
    }
    
    FitData <- BasicLineData
  }
  
  
  
  
  
  CombData <- FitData[[length(FitData)]]
  CombData$Report <- FitData[[1]]$Report
  
  for(i in 2:length(FitData)){
    CombData$Report <- rbind(CombData$Report, FitData[[i]]$Report)
  }
  
  Net <- list()
  TaxonList <- list()
  InfoData <- list()
  ProjPoints <- list()
  PCAPrGraph <- list()
  
  for(i in 1:length(FitData)){
    
    print(paste("Constructing accessory structures - round", i))
    
    Net[[i]] <- ConstructGraph(Results = FitData[[i]], DirectionMat = NULL, Thr = 0.05)
    TaxonList[[i]] <- getTaxonMap(Results = FitData[[i]], Data = Data)
    
    ProjPoints[[i]] <- projectPoints(Results = FitData[[i]], Data = Data, TaxonList = TaxonList[[i]],
                                     UseR = TRUE,
                                     method = 'PCALin', Dims = nDims, Debug = FALSE)
    
    if(i != length(FitData)){
      InfoData[[i]] <- plotPieNet(Results = FitData[[i]], Data = Data, NodeSizeMult = 4,
                                  Categories = Categories, PlotNet = FALSE,
                                  Graph = Net[[i]], TaxonList = TaxonList[[i]], LayOut =, Main = "Pincipal graph")
    } else {
      InfoData[[i]] <- plotPieNet(Results = FitData[[i]], Data = Data, NodeSizeMult = 4,
                                  Categories = Categories, PlotNet = TRUE,
                                  Graph = Net[[i]], TaxonList = TaxonList[[i]], LayOut =, Main = "Pincipal graph")
    }
    
    
    PCAPrGraph[[i]] <-  prcomp(FitData[[i]]$Nodes, retx = TRUE, center = FALSE, scale. = FALSE)
    
    # if(FitData[[i]]$Method == "CircleConfiguration"){
    #   RotatedData <- cbind(Data %*% PCAPrGraph[[i]]$rotation[,1:2], 1:nrow(Data) %in% NonG0Cell,
    #                        as.character(Categories))
    # } else {
    #   RotatedData <- cbind(Data %*% PCAPrGraph[[i]]$rotation[,1:2], rep(TRUE, nrow(Data)),
    #                        as.character(Categories))
    # }
    # 
    # colnames(RotatedData) <- c("PC1", "PC2", "NG0", "Cat")
    # 
    # RotatedData.DF <- data.frame(RotatedData)
    # RotatedData.DF$PC1 <- as.numeric(as.character(RotatedData.DF$PC1))
    # RotatedData.DF$PC2 <- as.numeric(as.character(RotatedData.DF$PC2))
    # RotatedData.DF$NG0 <- factor(RotatedData.DF$NG0, levels = c("TRUE", "FALSE"))
    # 
    # 
    # p <- ggplot(data.frame(RotatedData.DF), aes(x=PC1, y=PC2, alpha=NG0, colour=Cat)) + geom_point() +
    #   geom_point(data = data.frame(PCAPrGraph[[i]]$x[,1:2]), mapping = aes(x=PC1, y=PC2),
    #              inherit.aes = FALSE) +
    #   labs(title = paste("Round", i)) + scale_alpha_discrete("Fitted", range = c(1, .1))
    # 
    # for(j in 1:nrow(FitData[[i]]$Edges)){
    #   p <- p + geom_path(data = data.frame(PCAPrGraph[[i]]$x[FitData[[i]]$Edges[j,],1:2]),
    #                      mapping = aes(x = PC1, y = PC2), inherit.aes = FALSE)
    # }
    # 
    # print(p)
    
    if(PlotIntermediate & i != length(FitData)){
      if(FitData[[i]]$Method == "CircleConfiguration" & GraphType == 'Lasso'){
        ProjectOnPrincipalGraph(Nodes = FitData[[i]]$Nodes, Edges = FitData[[i]]$Edges, Points = Data,
                                UsedPoints = which(rownames(Data) %in% NonG0Cell), Categories = Categories, Title= paste("Round", i),
                                PCACenter = PCAProjCenter)
      } else {
        ProjectOnPrincipalGraph(Nodes = FitData[[i]]$Nodes, Edges = FitData[[i]]$Edges, Points = Data,
                                UsedPoints = NULL, Categories = Categories, Title= paste("Round", i),
                                PCACenter = PCAProjCenter)
      }
    } 
    
    if(i == length(FitData)){
      ProjectOnPrincipalGraph(Nodes = FitData[[i]]$Nodes, Edges = FitData[[i]]$Edges, Points = Data,
                              UsedPoints = NULL, Categories = Categories, Title= paste("Round", i),
                              PCACenter = PCAProjCenter)
    }
    
  }
  
  
  # legend("center", legend = levels(Categories), fill = InfoData$ColInfo[1:3], cex = 4)
  
  return(list(Data = Data, FiltExp = DataMat, Categories = Categories, PrinGraph = CombData,
              IntGrahs = FitData, Net = Net, TaxonList = TaxonList, PCAData = PCAData,
              nDims = nDims, ProjPoints = ProjPoints, PCAPrGraph = PCAPrGraph, NonG0Cell = NonG0Cell))
}

  

















################################################################################
#
# SelectGene ------------------------------------------------------
#
################################################################################

#' Filter Gene saccording to predefined properties
#'
#' @param BaseAnalysis a strure returned by Project and Compute
#' @param Mode scalar, a string indicating the algorithm to use. I can be
#' "VarPC" (Variance-based on genes the genes used to contruct the principal curve)
#' "PeakPC" (Peack-based on genes the genes used to contruct the principal curve)
#' "VarALL" (Variance-based on all genes)
#' "PeakALL" (Peack-based on all genes)
#' @param DistillThr numeric, a threshold used by the algorithm. It is a p-value for variance based selection
#' and a quantile for peack based selection
#' @param Topo string, a string describing the topology of the principal curve. It can be 'Lasso' or 'Circle'
#' @param FullExpression the matrix containing the full expression data (note that sample filtered during the construction of the principal graph will be removed)
#' @param LoesSpan the span parameter to use with the loes function
#' @param ExtMode 
#' @param IgnoreTail 
#' @param ZeroQuant 
#' @param Log 
#' @param MedianRatThr 
#' @param VarTestPVThr 
#'
#' @return
#' @export
#'
#' @examples
DistillGene <- function(BaseAnalysis, Mode = "VarPC", DistillThr = .2, ExtMode =1,
                        Topo = 'Lasso', IgnoreTail = FALSE, ZeroQuant = .75,
                        FullExpression = NULL, Log = FALSE, LoesSpan = .1,
                        MedianRatThr = .05, VarTestPVThr = 1e-5, FilterLow = FALSE) {
  
  if(!is.null(FullExpression) & Log){
    FullExpression <- log10(FullExpression + 1)
  } else {
    FullExpression <- FullExpression
  }
  
  
  # Get gene projection on the principal curve
  NodeOnGenes <- t(BaseAnalysis$PrinGraph$Nodes %*% t(BaseAnalysis$PCAData$rotation[,1:BaseAnalysis$nDims]))
  
  if(Mode == "VarPC" | Mode == "VarALL"){
    print("Selecting genes based on dispersion")
    
    # Find association of cells to nodes of the PG
    CellVertexAssociation <- rep(NA, nrow(BaseAnalysis$FiltExp))
    
    nNodes <- nrow(BaseAnalysis$PrinGraph$Nodes)
    
    for(i in 1:nNodes){
      CellVertexAssociation[ BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]][[i]] ] <- i
    }
    
    if(Mode == "VarALL"){
      print("Using variance over mean on node")
      
      # Only look at cells that are present in the principal graph structure
      FullExpression <- FullExpression[,rownames(BaseAnalysis$FiltExp)]
      
      # Split the data according to their assocaition to the different nodes
      SplitData <- lapply(split(1:length(CellVertexAssociation), CellVertexAssociation), function(x){FullExpression[,x]})
      
    } else {
      print("Using variance over the position of the principal curve. The absolute value ot the expression will be used")
      
      # Split the data according to their assocaition to the different nodes
      SplitData <- lapply(split(1:length(CellVertexAssociation), CellVertexAssociation), function(x){t(BaseAnalysis$FiltExp)[,x]})
      
    }

    # Only consider nodes with at least three cells projected
    ToUse <- paste(which(lapply(BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]], length) >= 3))
    
    # Consider only the circle if indicated
    if(Topo == 'Lasso' & IgnoreTail){
      TailPath <- GetLongestPath(Net = BaseAnalysis$Net[[length(BaseAnalysis$Net)]],
                                 Structure = 'Tail', Circular = FALSE)$VertNumb
      if(!is.null(TailPath)){
        ToUse <- setdiff(ToUse, TailPath[-length(TailPath)])
      }
    }
    
    SplitData <- SplitData[ToUse]
    NodeOnGenes <- NodeOnGenes[,as.numeric(ToUse)]
    
    if(Mode == "VarALL"){
      
      if(ExtMode == 1){
        
        print("Extended mode 1: median of the mean distances from the group means (0 means removed)")
        
        # Get the gene mean for each node
        MeanData <- lapply(SplitData, function(x){rowMeans(x)})
        
        # Get the difference between the data and the mean
        DiffData <- lapply(1:length(SplitData), function(i){apply( abs(SplitData[[i]] - MeanData[[i]]) , 1, mean) })
        
        # Put all together into matrices
        DistComb <- do.call(cbind, DiffData)
        MeanComb <- do.call(cbind, MeanData)
        
        tMat <- DistComb/MeanComb
        tMat[is.infinite(tMat)] <- NA
        
        StatData <- apply(tMat, 1, median, na.rm=TRUE)
        
      }
      
      if(ExtMode == 2){
        
        print("Extended mode 2: median coefficient of vartiation across groups (0 means removed)")
        
        # Get the gene mean for each node
        MeanData <- lapply(SplitData, function(x){rowMeans(x)})
        
        # Get the variance of the data
        VarData <- lapply(1:length(SplitData), function(i){apply(SplitData[[i]], 1, var)})
        
        # Put all together into matrices
        VarComb <- do.call(cbind, VarData)
        MeanComb <- do.call(cbind, MeanData)
        
        tMat <- VarComb/MeanComb
        tMat[is.infinite(tMat)] <- NA
        
        StatData <- apply(tMat, 1, median, na.rm=TRUE)
        
      }
      
      if(ExtMode == 3){
        
        print("Extended mode 3: mean IQR/median across groups (0 median removed)")
        
        IQRData <- lapply(SplitData, function(x){apply(x, 1, IQR)})
        MedianData <- lapply(SplitData, function(x){apply(x, 1, median)})
        
        IQRComb <- do.call(cbind, IQRData)
        MedianComb <- do.call(cbind, MedianData)
        
        IQRoMed <- IQRComb/MedianComb
        IQRoMed[!is.finite(IQRoMed)] <- NA
        
        StatData <- apply(IQRoMed, 1, mean, na.rm=TRUE)
        
      }
      
    }
    
    if(Mode == "VarPC"){
    
      if(ExtMode == 1){

        print("Extended mode 1: median of the mean distances from the principal curve (0 means removed)")
        
        # Get the difference between the data and the PC value
        DiffData <- lapply(1:length(SplitData), function(i){ apply(abs(SplitData[[i]] - NodeOnGenes[,i]), 1, mean) })

        # Put all together into matrices
        DistComb <- do.call(cbind, DiffData)

        # median normalized distance from the PC points
        tMat <- DistComb/abs(NodeOnGenes)
        tMat[is.infinite(tMat)] <- NA
        
        StatData <- apply(tMat, 1, median, na.rm=TRUE)

      }
      
      if(ExtMode == 2){

        print("Extended mode 2: median variance over the principal curve (0 variance removed)")
        
        # Get the variance of the data
        VarData <- lapply(1:length(SplitData), function(i){apply(SplitData[[i]], 1, var)})

        # Put all together into matrices
        VarComb <- do.call(cbind, VarData)

        # median coefficient of variation
        tMat <- VarComb/abs(NodeOnGenes)
        tMat[is.infinite(tMat)] <- NA
        
        StatData <- apply(tMat, 1, median, na.rm=TRUE)

      }
      
      if(ExtMode == 3){

        print("Extended mode 3: median IQR over the principal curve (0 IQR removed)")
        
        IQRData <- lapply(SplitData, function(x){apply(x, 1, IQR)})

        IQRComb <- do.call(cbind, IQRData)

        IQRoMed <- IQRComb/abs(NodeOnGenes)
        IQRoMed[!is.finite(IQRoMed)] <- NA

        StatData <- apply(IQRoMed, 1, median, na.rm=TRUE)

      }
      
    }
    
    KeepGenes <- names(which(StatData < DistillThr))
    KeepGenes <- KeepGenes[!is.na(KeepGenes)]
    
    print(paste(sum(StatData >= DistillThr, na.rm = TRUE),
                "genes filtered by the selected procedure"))
    print(paste(length(KeepGenes), "genes selected"))
    
  }
  
  if(Mode == "PeakPC" | Mode == "PeakALL"){
    print("Selecting genes based on the presence of single peaks")
    print("A random circular path will be chosen. Lasso tails will be ignored")

    CellVertexAssociation <- rep(NA, nrow(BaseAnalysis$FiltExp))

    nNodes <- nrow(BaseAnalysis$PrinGraph$Nodes)

    for(i in 1:nNodes){
      CellVertexAssociation[ BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]][[i]] ] <- i
    }

    AllPaths <- GetLongestPath(Net = BaseAnalysis$Net[[length(BaseAnalysis$Net)]],
                               Structure = Topo, Circular = TRUE)

    if(length(AllPaths$VertPath) > nrow(BaseAnalysis$PrinGraph$Nodes) + 1){
      PathToUse <- AllPaths$VertNumb[sample(1:nrow(AllPaths$VertNumb), 1),]
    } else {
      PathToUse <- AllPaths$VertNumb
    }

    if(Topo == 'Lasso'){
      TailPath <- GetLongestPath(Net = BaseAnalysis$Net[[length(BaseAnalysis$Net)]], Structure = 'Tail')$VertNumb
      if(!is.null(TailPath)){
        PathToUse <- PathToUse[!(PathToUse %in% TailPath[-length(TailPath)])]
      }
    }
    
    if(Mode == "PeakALL"){
      print("Using peak detection on all genes")
      
      # Only look at cells that are present in the principal graph structure
      FullExpression <- FullExpression[,rownames(BaseAnalysis$FiltExp)]
      
      # Split the data according to their assocaition to the different nodes
      SplitData <- lapply(split(1:length(CellVertexAssociation), CellVertexAssociation), function(x){FullExpression[,x]})
      
      ProjPoints <- BaseAnalysis$ProjPoints[[length(BaseAnalysis$ProjPoints)]]
      
      PathProjection <- OrderOnPath(PrinGraph = BaseAnalysis$PrinGraph, Path = PathToUse, PointProjections = ProjPoints)
      TotPathLen <- sum(PathProjection$PathLen)
      
    } else {
      print("Using peak detection on the principal curve")
      
      # Split the data according to their assocaition to the different nodes
      SplitData <- lapply(split(1:length(CellVertexAssociation), CellVertexAssociation), function(x){t(BaseAnalysis$FiltExp)[,x]})
      
    }
    
    if(Mode == "PeakPC"){
      
      if(ExtMode == 1){
        
        print("Extended mode 1: using the principal curve projections")
        
        # Remove the last point of the path (to prevent circularity)
        PathToUse <- PathToUse[-length(PathToUse)]
        
        RefExpr <- t(NodeOnGenes[,as.numeric(PathToUse)])
        
        BinMat <- RefExpr > apply(RefExpr, 1, quantile, DistillThr)
        
      }
      
      if(ExtMode == 2){
        
        print("Extended mode 2: using the means over the nodes")
        
        # Get the difference between the data and the PC value
        DiffData <- lapply(1:length(SplitData), function(i){ if(!is.null(dim(SplitData[[i]]))) rowMeans(SplitData[[i]]) else SplitData[[i]] })
        
        # Remove the last point of the path (to prevent circularity)
        PathToUse <- PathToUse[-length(PathToUse)]
        
        # Get the nodes without any point projected
        NonNodes <- which(!unlist(lapply(lapply(BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]], is.na), any)))
        FixedPathToUse <- PathToUse[PathToUse %in% paste(NonNodes)]
        
        # Put all together into matrices
        DistComb <- do.call(cbind, DiffData)
        colnames(DistComb) <- paste(NonNodes)
        DistComb <- DistComb[,FixedPathToUse]
        
        BinMat <- DistComb > apply(DistComb, 1, quantile, DistillThr)
        
      }
      
      print(paste(sum(rowSums(BinMat)==0 | rowSums(BinMat)==ncol(BinMat)), "genes filtered due to flatness"))
      BinMat <- BinMat[rowSums(BinMat)>0 & rowSums(BinMat)<ncol(BinMat),]
      
      Uni <- NULL
      Multi <- NULL

      for(i in 1:nrow(BinMat)){
        isUni <- FALSE
        if(all(BinMat[i,min(which(BinMat[i,])):max(which(BinMat[i,]))])){
          isUni <- TRUE
        }
        if(all(!BinMat[i,min(which(!BinMat[i,])):max(which(!BinMat[i,]))])){
          isUni <- TRUE
        }

        if(isUni){
          Uni <- c(Uni, i)
        } else {
          Multi <- c(Multi, i)
        }
      }
      
    }
    
    if(Mode == "PeakALL"){
      
      if(ExtMode == 1){
        
        print("Extended mode 1: using a smoother on ordered gene expression")
        
        print("Fitting loess smoothers. This can take some time")

        pb <- txtProgressBar(min = 0, max = nrow(FullExpression), style = 3)
        
        PathProjection$PositionOnPath == 0

        X <- c(PathProjection$PositionOnPath - TotPathLen,
               PathProjection$PositionOnPath,
               PathProjection$PositionOnPath + TotPathLen)
        
        Y <- cbind(FullExpression[,order(PathProjection$PositionOnPath)],
                   FullExpression[,order(PathProjection$PositionOnPath)],
                   FullExpression[,order(PathProjection$PositionOnPath)])
        
        AllSmooth <- sapply(X = as.list(1:nrow(FullExpression)), function(i){
          setTxtProgressBar(pb, value = i)
          predict(loess(unlist(Y[i,]) ~ X, span = LoesSpan), newdata = data.frame(X = cumsum(PathProjection$PathLen)))
        })
        AllSmooth <- t(AllSmooth)
        
        BinMat <- AllSmooth >= apply(AllSmooth, 1, quantile, DistillThr)
        
      }
      
      if(ExtMode == 2){
        
        print("Extended mode 2: using the means over the nodes")
        
        # Get the difference between the data and the PC value
        DiffData <- lapply(1:length(SplitData), function(i){ if(!is.null(dim(SplitData[[i]]))) rowMeans(SplitData[[i]]) else SplitData[[i]] })
        
        # Remove the last point of the path (to prevent circularity)
        PathToUse <- PathToUse[-length(PathToUse)]
        
        # Get the nodes without any point projected
        NonNodes <- which(!unlist(lapply(lapply(BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]], is.na), any)))
        FixedPathToUse <- PathToUse[PathToUse %in% paste(NonNodes)]
        
        # Put all together into matrices
        DistComb <- do.call(cbind, DiffData)
        colnames(DistComb) <- paste(NonNodes)
        DistComb <- DistComb[,FixedPathToUse]
        
        BinMat <- DistComb > apply(DistComb, 1, quantile, DistillThr)
        
      }
      
      print(paste(sum(rowSums(BinMat)==1 | rowSums(BinMat)==ncol(BinMat)), "genes filtered due to flatness"))
      BinMat <- BinMat[rowSums(BinMat)>1 & rowSums(BinMat)<ncol(BinMat),]
      
      Uni <- NULL
      Multi <- NULL
      
      for(i in 1:nrow(BinMat)){
        isUni <- FALSE
        if(all(BinMat[i,min(which(BinMat[i,])):max(which(BinMat[i,]))])){
          isUni <- TRUE
        }
        if(all(!BinMat[i,min(which(!BinMat[i,])):max(which(!BinMat[i,]))])){
          isUni <- TRUE
        }
        
        if(isUni){
          Uni <- c(Uni, i)
        } else {
          Multi <- c(Multi, i)
        }
      }
      
    }
     
    KeepGenes <- rownames(BinMat)[Uni]
    
    print(paste(length(Multi), "genes filtered due to multimodality"))
    print(paste(length(KeepGenes), "genes selected"))
    
     
  }
  
  if(FilterLow){
    # Get the median expresion value 
    MedianPCExp <- apply(NodeOnGenes, 1, quantile, ZeroQuant)
    
    Jumps <- (sort(MedianPCExp)[-1] - sort(MedianPCExp)[-length(MedianPCExp)])
    BigJmps <- which((Jumps-mean(Jumps))/sd(Jumps)>3)
    
    PCnOK <- names(Jumps)[1:(BigJmps[1]-1)]
    PCOK <- names(Jumps)[(BigJmps[1]):length(Jumps)]
    
    if(length(PCOK) > 10 & length(PCnOK) > 10){
      MedianRat <- median(MedianPCExp[PCnOK])/median(MedianPCExp)
      VarPopPV <- var.test(MedianPCExp[PCnOK], MedianPCExp[PCOK], alternative = "less")$p.value
      
      if(MedianRat < MedianRatThr & VarPopPV < VarTestPVThr){
        return(intersect(KeepGenes, PCOK))
      }
    }
  }

  
  return(KeepGenes)
    
}


####### Extra code



# CombinedMat <- cbind(DiffExpData, scale(BaseAnalysis$FiltExp, center = TRUE))
# 
# VarDiff <- apply(CombinedMat, 2, function(x) {wilcox.test(x[1:(length(x)/2)], x[((length(x)/2)+1):length(x)],
#                                                        alternative = "less")$p.value} )
# KeepGenes <- names(which(VarDiff<DistillThr))

# hist( abs(apply(DiffExpData, 2, median)) / apply(MatrixExpansion, 2, median) )



# par(mfcol=c(3,3))
# 
# for(i in 1:length(KeepGenes)) boxplot(list(DiffExpData[, KeepGenes[i]],
#                                            scale(BaseAnalysis$FiltExp[, KeepGenes[i]], center = TRUE)))



# if(Mode == "VarALL"){
#   print("Selecting genes with the smallest fluctuations")
#   
#   # AllPaths <- GetLongestPath(Net = BaseAnalysis$Net[[length(BaseAnalysis$Net)]],
#   #                            Structure = Topo, Circular = TRUE)
#   # 
#   # ProjPoints <- BaseAnalysis$ProjPoints[[length(BaseAnalysis$ProjPoints)]]
#   # 
#   # if(length(AllPaths$VertPath) > nrow(BaseAnalysis$PrinGraph$Nodes) + 1){
#   #   PathToUse <- AllPaths$VertNumb[sample(1:nrow(AllPaths$VertNumb), 1),]
#   # } else {
#   #   PathToUse <- AllPaths$VertNumb
#   # }
#   # 
#   # PathProjection <- OrderOnPath(PrinGraph = BaseAnalysis$PrinGraph, Path = PathToUse, PointProjections = ProjPoints)
#   # TotPathLen <- sum(PathProjection$PathLen)
#   
#   
#   # Find association of cells to nodes of the PG
#   CellVertexAssociation <- rep(NA, nrow(BaseAnalysis$FiltExp))
#   
#   nNodes <- nrow(BaseAnalysis$PrinGraph$Nodes)
#   
#   for(i in 1:nNodes){
#     CellVertexAssociation[ BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]][[i]] ] <- i
#   }
#   
#   # Only look at cells that are present in the principal graph structure
#   FullExpression <- FullExpression[,rownames(BaseAnalysis$FiltExp)]
#   
#   # Split the data according to their assocaition to the different nodes
#   SplitData <- lapply(split(1:length(CellVertexAssociation), CellVertexAssociation), function(x){FullExpression[,x]})
#   
#   # Only consider nodes with at least three cells projected
#   ToUse <- paste(which(lapply(BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]], length) >= 3))
#   SplitData <- SplitData[ToUse]
#   
#   if(ExtMode == 1){
#     
#     # Get the gene mean for each node
#     MeanData <- lapply(SplitData, function(x){rowMeans(x)})
#     
#     # Get the difference between the data and the mean
#     DiffData <- lapply(1:length(SplitData), function(i){SplitData[[i]] - MeanData[[i]]})
#     
#     # Put all together into matrices
#     DistComb <- do.call(cbind, DiffData)
#     MeanComb <- do.call(cbind, MeanData)
#     
#     StatData <- apply(abs(DistComb)/MeanComb, 1, median)
#     
#   }
#   
#   
#   if(ExtMode == 2){
#     
#     # Get the gene mean for each node
#     MeanData <- lapply(SplitData, function(x){rowMeans(x)})
#     
#     # Get the variance of the data
#     VarData <- lapply(1:length(SplitData), function(i){apply(SplitData[[i]], 1, var)})
#     
#     # Put all together into matrices
#     VarComb <- do.call(cbind, VarData)
#     MeanComb <- do.call(cbind, MeanData)
#     
#     StatData <- apply(VarComb/MeanComb, 1, median)
#     
#   }
#   
#   
#   if(ExtMode == 3){
#     
#     IQRData <- lapply(SplitData, function(x){if(is.null(dim(x))) rep(0, length(x)) else apply(x, 1, IQR)})
#     MedianData <- lapply(SplitData, function(x){if(is.null(dim(x))) x else apply(x, 1, median)})
#     
#     IQRComb <- do.call(cbind, IQRData)
#     MedianComb <- do.call(cbind, MedianData)
#     
#     IQRoMed <- IQRComb/MedianComb
#     IQRoMed[!is.finite(IQRoMed)] <- NA
#     
#     StatData <- apply(IQRoMed, 1, mean, na.rm=TRUE)
#     
#   }
#   
#   
#   
#   # print("Fitting loess smoothers")
#   # 
#   # pb <- txtProgressBar(min = 0, max = nrow(FullExpression), style = 3)
#   # 
#   # AllSmooth <- lapply(X = as.list(1:nrow(FullExpression)), function(i){
#   #   setTxtProgressBar(pb, value = i)
#   #   return(loess(unlist(FullExpression[i,]) ~ PathProjection$PositionOnPath, span = LoesSpan))
#   # })
#   #  
#   # AllResiduals <- lapply(lapply(AllSmooth, "[[", "residuals"), var)
#   # 
#   # print("Testing the distribution of the residuals")
#   # 
#   # pb1 <- txtProgressBar(min = 0, max = nrow(FullExpression), style = 3)
#   # 
#   # AllWT <- lapply(as.list(1:nrow(FullExpression)), function(i){
#   #   setTxtProgressBar(pb1, value = i)
#   #   wilcox.test(unlist(FullExpression[i,]) - mean(unlist(FullExpression[i,])), AllSmooth[[i]]$residuals, alternative = "greater" )
#   # })
#   
#   # AllMedRat <- lapply(as.list(1:nrow(FullExpression)), function(i){
#   #   setTxtProgressBar(pb1, value = i)
#   #   median(AllSmooth[[i]]$residuals) - median(unlist(FullExpression[i,]) - mean(unlist(FullExpression[i,])))
#   # })
#   
#   # AllWT <- AllPV
#   
#   # PVVect <- unlist(lapply(AllWT, "[[", "p.value"))
#   # names(PVVect) <- rownames(FullExpression)
#   
#   print(paste(sum(StatData > DistillThr, na.rm = TRUE), "genes excluded due to thresholding"))
#   
#   KeepGenes <- rownames(FullExpression)[StatData < DistillThr]
#   KeepGenes <- KeepGenes[!is.na(KeepGenes)]
#   
#   return(KeepGenes)
#   
#   # i <- 10
#   # 
#   # boxplot(list(
#   #   (unlist(FullExpression[order(PVVect, decreasing = FALSE)[i],]) - mean(unlist(FullExpression[order(PVVect, decreasing = FALSE)[i],]))),
#   #   AllSmooth[[order(PVVect, decreasing = FALSE)[i]]]$residuals
#   # ), main = sort(PVVect)[i])
#   
# }




















# if(Mode == "PeakPC"){
#   print("Selecting genes with more defined single peaks on the principal curve")
#   
#   CellVertexAssociation <- rep(NA, nrow(BaseAnalysis$FiltExp))
#   
#   nNodes <- nrow(BaseAnalysis$PrinGraph$Nodes)
#   
#   for(i in 1:nNodes){
#     CellVertexAssociation[ BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]][[i]] ] <- i
#   }
#   
#   NodeOnGenes <- BaseAnalysis$PrinGraph$Nodes %*% t(BaseAnalysis$PCAData$rotation[,1:BaseAnalysis$nDims])
#   
#   AllPaths <- GetLongestPath(Net = BaseAnalysis$Net[[length(BaseAnalysis$Net)]],
#                              Structure = Topo, Circular = TRUE)
#   
#   if(length(AllPaths$VertPath) > nrow(BaseAnalysis$PrinGraph$Nodes) + 1){
#     PathToUse <- AllPaths$VertNumb[sample(1:nrow(AllPaths$VertNumb), 1),]
#   } else {
#     PathToUse <- AllPaths$VertNumb
#   }
#   
#   if(Topo == 'Lasso'){
#     Tail <- GetLongestPath(Net = BaseAnalysis$Net[[length(BaseAnalysis$Net)]], Structure = 'tail')$VertNumb
#     PathToUse <- PathToUse[!(PathToUse %in% Tail[-length(Tail)])]
#   }
#   
#   NodeOnGenes <- NodeOnGenes[as.integer(PathToUse),]
#   
#   BinMat <- t(NodeOnGenes) >= apply(NodeOnGenes, 2, quantile, DistillThr)
#   
#   BinMat <- BinMat[rowSums(BinMat)>0, ]
#   
#   Uni <- NULL
#   Multi <- NULL
#   
#   for(i in 1:nrow(BinMat)){
#     isUni <- FALSE
#     if(all(BinMat[i,min(which(BinMat[i,])):max(which(BinMat[i,]))])){
#       isUni <- TRUE
#     }
#     if(all(!BinMat[i,min(which(!BinMat[i,])):max(which(!BinMat[i,]))])){
#       isUni <- TRUE
#     }
#     
#     if(isUni){
#       Uni <- c(Uni, i)
#     } else {
#       Multi <- c(Multi, i)
#     }
#   }
#   
#   # gplots::heatmap.2(1*(BinMat), Colv = FALSE, dendrogram = 'row')
#   
#   KeepGenes <- rownames(BinMat)[Uni]
#   
#   return(KeepGenes)
#   
# }



# if(Mode == "PeakALL"){
#   print("Selecting genes with more defined single peaks in the smoother")
#   
#   AllPaths <- GetLongestPath(Net = BaseAnalysis$Net[[length(BaseAnalysis$Net)]],
#                              Structure = Topo, Circular = TRUE)
#   
#   ProjPoints <- BaseAnalysis$ProjPoints[[length(BaseAnalysis$ProjPoints)]]
#   
#   if(length(AllPaths$VertPath) > nrow(BaseAnalysis$PrinGraph$Nodes) + 1){
#     PathToUse <- AllPaths$VertNumb[sample(1:nrow(AllPaths$VertNumb), 1),]
#   } else {
#     PathToUse <- AllPaths$VertNumb
#   }
#   
#   if(Topo == 'Lasso'){
#     Tail <- GetLongestPath(Net = BaseAnalysis$Net[[length(BaseAnalysis$Net)]], Structure = 'tail')$VertNumb
#     PathToUse <- PathToUse[!(PathToUse %in% Tail[-length(Tail)])]
#   }
#   
#   # PathProjection <- OrderOnPath(PrinGraph = BaseAnalysis$PrinGraph, Path = PathToUse, PointProjections = ProjPoints)
#   # TotPathLen <- sum(PathProjection$PathLen)
#   
#   FullExpression <- FullExpression[,rownames(BaseAnalysis$FiltExp)]
#   
#   CellVertexAssociation <- rep(NA, nrow(BaseAnalysis$FiltExp))
#   
#   nNodes <- nrow(BaseAnalysis$PrinGraph$Nodes)
#   
#   for(i in 1:nNodes){
#     CellVertexAssociation[ BaseAnalysis$TaxonList[[length(BaseAnalysis$TaxonList)]][[i]] ] <- i
#   }
#   
#   SplitData <- lapply(split(1:length(CellVertexAssociation), CellVertexAssociation), function(x){FullExpression[,x]})
#   MeanData <- lapply(SplitData, function(x){if(is.null(dim(x))) x else rowMeans(x)})
#   MeanComb <- do.call(cbind, MeanData)
#   
#   ReordedMeans <- MeanComb[,PathToUse[PathToUse %in% colnames(MeanComb)]]
#   
#   BinMat <- ReordedMeans >= apply(ReordedMeans, 1, quantile, DistillThr)
#   
#   # print("Fitting loess smoothers")
#   # 
#   # pb <- txtProgressBar(min = 0, max = nrow(FullExpression), style = 3)
#   # 
#   # X <- PathProjection$PositionOnPath
#   # X[X == 0 | X == TotPathLen] <- 0
#   # OutPath <- is.na(PathProjection$PositionOnPath)
#   # X <- X[!OutPath]
#   # X <- c(X-TotPathLen, X, X+TotPathLen)
#   # 
#   # AllSmooth <- lapply(X = as.list(1:nrow(FullExpression)), function(i){
#   #   setTxtProgressBar(pb, value = i)
#   # 
#   #   Y <- unlist(FullExpression[i,!OutPath])
#   #   Y <- c(Y, Y, Y)
#   #   
#   #   return(loess(Y ~ X, span = LoesSpan))
#   # })
#   # 
#   # print("Obtaining smoothed values")
#   # 
#   # pb1 <- txtProgressBar(min = 0, max = nrow(FullExpression), style = 3)
#   # 
#   # AllPos <- cumsum(PathProjection$PathLen)
#   # 
#   # SmootherData <- sapply(as.list(1:length(AllSmooth)), function(i) {
#   #   setTxtProgressBar(pb1, value = i)
#   #   predict(AllSmooth[[i]], newdata = data.frame(X=AllPos))
#   # })
#   # 
#   # colnames(SmootherData) <- rownames(FullExpression)
# 
#   # BinMat <- t(SmootherData) >= apply(SmootherData, 2, quantile, DistillThr)
#   
#   BinMat <- BinMat[rowSums(BinMat)>0 & rowSums(BinMat)<ncol(BinMat), ]
#   
#   Uni <- NULL
#   Multi <- NULL
#   UniBroad <- NULL
#   
#   for(i in 1:nrow(BinMat)){
#     isUni <- FALSE
#     if(all(BinMat[i,min(which(BinMat[i,])):max(which(BinMat[i,]))])){
#       isUni <- TRUE
#     }
#     if(all(!BinMat[i,min(which(!BinMat[i,])):max(which(!BinMat[i,]))])){
#       isUni <- TRUE
#     }
#     
#     if(isUni){
#       Uni <- c(Uni, i)
#       UniBroad <- c(UniBroad, sum(BinMat[i,]))
#     } else {
#       Multi <- c(Multi, i)
#     }
#   }
#   
#   # gplots::heatmap.2(1*(BinMat), Colv = FALSE, dendrogram = 'row')
#   
#   KeepGenes <- rownames(BinMat)[Uni[UniBroad <= RatioPeaks*ncol(BinMat)]]
#   
#   return(KeepGenes)
#   
# }












################################################################################
#
# Converge on genes ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param ExpData 
#' @param StartGeneSet 
#' @param Mode 
#' @param DistillThr 
#' @param ExtMode 
#' @param StopCrit 
#' @param OutThr 
#' @param nNodes 
#' @param VarThr 
#' @param Categories 
#' @param GraphType 
#' @param PlanVarLimit 
#' @param PlanVarLimitIC 
#' @param ForceLasso 
#' @param LassoCircInit 
#' @param MinBranDiff 
#' @param Log 
#' @param Filter 
#' @param MinProlCells 
#' @param DipPVThr 
#' @param PCACenter 
#' @param PCAProjCenter 
#' @param PlotIntermediate 
#'
#' @return
#' @export
#'
#' @examples
ConvergeOnGenes <- function(ExpData, StartGeneSet, Mode = "VarALL", DistillThr = .5, ExtMode = 1, StopCrit = .99, MaxRounds = 25, FastReduce = FALSE,
                            OutThr = 5, nNodes = 40, VarThr = .99, Categories = NULL, UpdateG0Cells = TRUE, PCAFilter = TRUE, OutThrPCA = 5,
                            GraphType = 'Circle', IgnoreTail = FALSE, PlanVarLimit = .90, PlanVarLimitIC = .95, ForceLasso = FALSE,
                            InitStructNodes = 20, MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 20,
                            DipPVThr = 1e-4, PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE,
                            StartQuant = .5, QuantInc = .01, EstProlif = "MeanPerc", PlotDebug = FALSE) {
  
  StartStruct <- ProjectAndCompute(DataSet = ExpData, GeneSet = StartGeneSet, OutThr = OutThr, nNodes = nNodes,
                                         VarThr = VarThr, Categories = Categories, GraphType = GraphType, PlanVarLimit = PlanVarLimit,
                                         PlanVarLimitIC = PlanVarLimitIC, ForceLasso = ForceLasso, InitStructNodes = InitStructNodes,
                                         MinBranDiff = MinBranDiff, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                         PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = PlotIntermediate,
                                         QuaThr = StartQuant, EstProlif = EstProlif, PCAFilter = PCAFilter, OutThrPCA = OutThrPCA,
                                   PlotDebug = PlotDebug)
  ConsideredStruct <- StartStruct
  
  FilteredGenes <- StartGeneSet
  Converged <- FALSE
  RoundCount <- 0
  GeneNumber <- NULL
  
  while(!Converged){
    
    RoundCount <- RoundCount + 1
    
    OldFiltered <- FilteredGenes
    
    FilteredGenes <- DistillGene(BaseAnalysis = ConsideredStruct,
                                 Mode = Mode, DistillThr = DistillThr, ExtMode = ExtMode,
                                 Topo = GraphType, IgnoreTail = IgnoreTail, FullExpression = ExpData, Log = Log)
    
    if(is.null(FilteredGenes)){
      print("No gene found")
      return(NULL)
    }
    
    if(length(FilteredGenes) == 0){
      print("No gene found")
      return(NULL)
    }
    
    GeneNumber <- c(GeneNumber, length(FilteredGenes))
    
    print(paste(length(FilteredGenes), "genes selected"))
    
    if(length(intersect(FilteredGenes, OldFiltered))/length(OldFiltered) > StopCrit){
      plot(GeneNumber, xlab="Iterations", ylab = "Number of genes")
      return(FilteredGenes)
    }
    
    if(RoundCount > MaxRounds){
      print("Max number of iterations reached, returning")
      plot(GeneNumber, xlab="Iterations", ylab = "Number of genes")
      return(FilteredGenes)
    }
    
    if(FastReduce){
      ExpData <- ExpData[rownames(ExpData) %in% FilteredGenes,]
    }
    
    if(UpdateG0Cells){
      ConsideredStruct <- ProjectAndCompute(DataSet = ExpData, GeneSet = FilteredGenes, OutThr = OutThr, nNodes = nNodes,
                                       VarThr = VarThr, Categories = Categories, GraphType = GraphType, PlanVarLimit = PlanVarLimit,
                                       PlanVarLimitIC = PlanVarLimitIC, ForceLasso = ForceLasso, InitStructNodes = InitStructNodes,
                                       MinBranDiff = MinBranDiff, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                       PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = PlotIntermediate,
                                       QuaThr = StartQuant + RoundCount*QuantInc, NonG0Cell = NULL, EstProlif = EstProlif,
                                       PCAFilter = PCAFilter, OutThrPCA = OutThrPCA, PlotDebug = PlotDebug)
    } else {
      ConsideredStruct <- ProjectAndCompute(DataSet = ExpData, GeneSet = FilteredGenes, OutThr = OutThr, nNodes = nNodes,
                                       VarThr = VarThr, Categories = Categories, GraphType = GraphType, PlanVarLimit = PlanVarLimit,
                                       PlanVarLimitIC = PlanVarLimitIC, ForceLasso = ForceLasso, InitStructNodes = InitStructNodes,
                                       MinBranDiff = MinBranDiff, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                       PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = PlotIntermediate,
                                       QuaThr = StartQuant + RoundCount*QuantInc, NonG0Cell = StartStruct$NonG0Cell, EstProlif = EstProlif,
                                       PCAFilter = PCAFilter, OutThrPCA = OutThrPCA, PlotDebug = PlotDebug)
    }

    
  }
  
}









# Extra code --------------------------------------------------------------

# if(EstProlif =="DipQuant"){
#   
#   print("Finding most significant multimodality thr")
#   
#   CheckVals <- seq(from = .25, by=.025, to = .75)
#   DipPV <- NULL
#   
#   for(i in CheckVals){
#     print(paste("Looking for dip at", i))
#     DipPV <- c(DipPV, diptest::dip.test(apply(RankedData, 2, quantile, i),
#                                         simulate.p.value = TRUE, B = 5000)$p.value)
#     
#     if(DipPV[length(DipPV)] < DipPVThr){
#       if(PlotDebug){
#         hist(apply(RankedData, 2, quantile, i),
#              main = "Gene expression ranking",
#              xlab = paste("q", i, "of gene expression by cell"))
#       }
#       
#     }
#   }
#   
#   IDxs <- which(DipPV < DipPVThr)
#   Association <- NULL
#   
#   for(i in IDxs){
#     
#     print(paste("Using mixed gaussian model on q", CheckVals[i]))
#     
#     MixedModel <- try(mixtools::normalmixEM(apply(RankedData, 2, quantile, CheckVals[i]),
#                                             maxit = 1000, maxrestarts = 100), TRUE)
#     
#     if(!is.list(MixedModel)){
#       next()
#     }
#     
#     if(ncol(MixedModel$posterior) > 2){
#       print("Multiple populations detected, moving on ...")
#       next()
#     }
#     
#     if(ncol(MixedModel$posterior) == 1){
#       print("Only one population detected, moving on ...")
#       next()
#     }
#     
#     print("Looking for quantile separation")
#     
#     QA <- quantile(MixedModel$x[MixedModel$posterior[,1]>.5], probs = c(.25, .75))
#     QB <- quantile(MixedModel$x[MixedModel$posterior[,2]>.5], probs = c(.25, .75))
#     
#     if(any(is.na(QA)) | any(is.na(QB))){
#       print("Too few cells in at least one population, moving on ...")
#       next()
#     }
#     
#     if( (min(QA) < min(QB) & max(QA) < min(QB)) |
#         (min(QB) < min(QA) & max(QB) < min(QA)) ){
#       
#       print("Quantile separation detected")
#       
#       Association <- rbind(Association, MixedModel$posterior[,which.max(MixedModel$mu)])
#       
#       if(PlotDebug){
#         boxplot(MixedModel$x[MixedModel$posterior[,1]>.5], MixedModel$x[MixedModel$posterior[,2]>.5],
#                 main=paste(CheckVals[i], "Quantile separated"))
#       }
#     } else {
#       if(PlotDebug){
#         boxplot(MixedModel$x[MixedModel$posterior[,1]>.5], MixedModel$x[MixedModel$posterior[,2]>.5],
#                 main=paste(CheckVals[i], "Quantile non-separated"))
#       }
#       
#     }
#     
#   }
#   
#   
#   if(length(Association)>ncol(RankedData)){
#     NonG0Cell <-  which(apply(Association, 2, min) > .5)
#     G0Cell <-  which(apply(Association, 2, min) <= .5)
#   } else {
#     if(!is.null(Association)){
#       NonG0Cell <- which(Association > .5)
#       G0Cell <- which(Association <= .5)
#     } else {
#       NonG0Cell <- NULL
#       G0Cell <- NULL
#     }
#     
#   }
#   
# }
# 










################################################################################
#
# PlotOnStages ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param Structure 
#' @param TaxonList 
#' @param Categories 
#' @param PrinGraph 
#' @param Net 
#' @param SelThr 
#' @param ComputeOverlaps 
#' @param RotatioMatrix 
#' @param ExpData 
#'
#' @return
#' @export
#'
#' @examples
PlotOnStages <- function(Structure, TaxonList, Categories, PrinGraph, Net, SelThr = NULL,
                         ComputeOverlaps = TRUE, RotatioMatrix, PointProjections, ExpData, nGenes = 10) {
  
  if(!is.factor(Categories)){
    stop("Categories must be a factor")
  }
  
  if(is.null(SelThr)){
    SelThr <- 1.1*(1/length(levels(Categories)))
  }
  
  # Structure <- "Circle"
  # TaxonList <- BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]]
  # Categories <- BuettInfo$FinalStruct$Categories
  # Net <- BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]]
  # RotatioMatrix <- BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims]
  # PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]]
  # ExpData <- BuettInfo$FinalStruct$FiltExp
  
  Nodes <- PrinGraph$Nodes
  nNodes <- nrow(PrinGraph$Nodes)
  
  # Find the best Staging if circular 
  
  TaxVect <- rep(NA, nNodes)
  
  for(i in 1:length(TaxonList)){
    TaxVect[TaxonList[[i]]] <- i 
  }
  
  TaxVect <- factor(TaxVect, levels = paste(1:nNodes))
  
  TB <- table(Categories, TaxVect)
  # colnames(TB)
  
  print("Step I - Finding the best path")
  
  StepIDone <- FALSE
  
  if(Structure == 'Circle'){
    
    StepIDone <- TRUE
    
    print("Getting all circular subisomorphisms")
    
    AllPaths <- GetLongestPath(Net = Net, Structure = Structure, Circular = TRUE)
    
    SummInfo <- NULL
    
    for(i in 1:nrow(AllPaths$VertNumb)){
      SelPath <- AllPaths$VertNumb[i,]
      # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
      SelPath <- SelPath[-length(SelPath)]
      
      Reordered <- TaxVect
      for(j in 1:length(SelPath)){
        Reordered[as.character(TaxVect) == SelPath[j]] <- j
      }
      
      NumReord <- as.numeric(as.character(Reordered))
      
      AGG <- aggregate(NumReord, by = list(Categories), median)
      AGG2 <- aggregate(NumReord, by = list(Categories), min)
      AGG3 <- aggregate(NumReord, by = list(Categories), max)
      
      Sorted <- FALSE
      
      if(S4Vectors::isSorted(AGG[,2])){
        SummInfo <- rbind(SummInfo,
                          c(i, 1, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                            AGG2[1,2], AGG[1,2])
        )
        Sorted <- TRUE
        # boxplot(NumReord ~ Categories, main = i)
      }
      
      if(S4Vectors::isSorted(rev(AGG[,2]))){
        SummInfo <- rbind(SummInfo,
                          c(i, 2, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                            nNodes - AGG3[1,2] + 1,
                            nNodes - AGG[1,2] + 1)
        )
        Sorted <- TRUE
        # boxplot(NumReord ~ Categories, main = paste(i, "rev"))
      }
      
      if(!Sorted){
        
        if(cor(AGG[,2], 1:nrow(AGG))>0){
          SummInfo <- rbind(SummInfo,
                            c(i, 3, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                              AGG2[1,2], AGG[1,2])
          )
        } else {
          SummInfo <- rbind(SummInfo,
                            c(i, 4, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                              nNodes - AGG3[1,2] + 1,
                              nNodes - AGG[1,2] + 1)
          )
        }
        
      }
    }
    
    print(SummInfo)
    
    print("Finding the best path")
    
    if(any(SummInfo[,2] %in% 1:2)){
      SummInfo <- SummInfo[SummInfo[,2] %in% 1:2, ]
    } else {
      print("Unable to find strongly consecutive stages")
    }
    
    Selected <- SummInfo[SummInfo[, 3] == min(SummInfo[, 3]), ]
    
    print(Selected)
    
    if(length(Selected)>5){
      Selected <- Selected[which.min(Selected[,4]),]
    }
    
    print(Selected)
    
    SelPath <- AllPaths$VertNumb[Selected[1],]
    # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
    SelPath <- SelPath[-length(SelPath)]
    
    if(Selected[2] %in% c(2, 4)){
      SelPath <- rev(SelPath)
    }
    
    Reordered <- TaxVect
    for(j in 1:length(SelPath)){
      Reordered[as.character(TaxVect) == SelPath[j]] <- j
    }
    
    NumReord <- as.numeric(as.character(Reordered))
    
    boxplot(NumReord ~ Categories)
    
  }
  
  if(Structure == 'Line'){
    
    StepIDone <- TRUE
    
    print("Getting all line subisomorphisms")
    
    AllPaths <- GetLongestPath(Net = Net, Structure = Structure, Circular = FALSE)
    
    SummInfo <- NULL
    
    for(i in 1:nrow(AllPaths$VertNumb)){
      SelPath <- AllPaths$VertNumb[i,]
      # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
      # SelPath <- SelPath[-length(SelPath)]
      
      Reordered <- TaxVect
      for(j in 1:length(SelPath)){
        Reordered[TaxVect == SelPath[j]] <- j
      }
      
      NumReord <- as.numeric(as.character(Reordered))
      
      AGG <- aggregate(NumReord, by = list(Categories), median)
      AGG2 <- aggregate(NumReord, by = list(Categories), min)
      AGG3 <- aggregate(NumReord, by = list(Categories), max)
      
      if(S4Vectors::isSorted(AGG[,2])){
        SummInfo <- rbind(SummInfo,
                          c(i, 1, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                            AGG2[1,2], AGG[1,2])
        )
        
        # boxplot(NumReord ~ Categories, main = i)
        
      }
      
      
      if(S4Vectors::isSorted(rev(AGG[,2]))){
        SummInfo <- rbind(SummInfo,
                          c(i, 2, summary(aov(NumReord ~ Categories))[[1]][1,"Pr(>F)"],
                            nNodes - AGG3[1,2] + 1,
                            nNodes - AGG[1,2] + 1)
        )
        # boxplot(Reordered ~ Categories, main = paste(i, "rev"))
      }
    }
    
    print("Finding the best path")
    
    Selected <- SummInfo[SummInfo[, 3] == min(SummInfo[, 3]), ]
    
    if(length(Selected)>5){
      Selected <- Selected[which.min(Selected[,4]),]
    }
    
    SelPath <- AllPaths$VertNumb[Selected[1],]
    # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
    SelPath <- SelPath[-length(SelPath)]
    
    if(Selected[2] == 2){
      SelPath <- rev(SelPath)
    }
    
    Reordered <- TaxVect
    for(j in 1:length(SelPath)){
      Reordered[as.character(TaxVect) == SelPath[j]] <- j
    }
    
    NumReord <- as.numeric(as.character(Reordered))
    
    boxplot(NumReord ~ Categories)
    
  }
  
  if(!StepIDone){
    stop("Unsupported Structure")
  }
  
  
  print("Step II")
  
  Empty <- which(unlist(lapply(lapply(TaxonList, is.na), any)))
  
  ExtPath <- SelPath

  if(Structure == 'Circle'){
    ExtPath <- c(ExtPath, ExtPath[1])
  }
  
  ExtendedTB <- TB[,ExtPath]
  
  barplot(t(t(ExtendedTB)/colSums(ExtendedTB)), col = rainbow(nrow(TB)), beside = TRUE, las = 2,
          legend.text = rownames(ExtendedTB), args.legend = list(x = "top", fill=rainbow(nrow(TB))),
          ylim = c(0, 1.25), yaxt = "n")
  axis(2, seq(from=0, to=1, by=.25), las=2)
  
  # SelPath <- SelPath[-length(SelPath)]
  SelPathSTG <- rep(NA, length(SelPath))
  
  
  PercMat <- t(t(TB)/colSums(TB))
  PercMat[is.na(PercMat)] <- 0
  PercMat[,colSums(TB)<3] <- 0
  PercMat <- PercMat[,ExtPath]
  BinPercMat <- (PercMat > SelThr)
  
  # if(!any(BinPercMat[,1])){
    # BinPercMat[,1] <- BinPercMat[,2]
  # }
  
  if(Structure == 'Circle'){
    
    for(j in 1:nrow(BinPercMat)){
      if(BinPercMat[j,2] & BinPercMat[j,ncol(BinPercMat)]){
        BinPercMat[j,1] <- TRUE
      }
    }
    
  }
  
  
  for(i in 2:(ncol(BinPercMat)-1)){
    for(j in 1:nrow(BinPercMat)){
      if(BinPercMat[j,i-1] & BinPercMat[j,i+1]){
        BinPercMat[j,i] <- TRUE
      }
    }
  }
  
  
  
  if(Structure == 'Circle'){
    
    for(j in 1:nrow(BinPercMat)){
      if(BinPercMat[j,1] & BinPercMat[j,ncol(BinPercMat)-1]){
        BinPercMat[j,ncol(BinPercMat)] <- TRUE
      }
    }
    
  }
  
  
  
  
  
  
  
  
  
  
  
  if(Structure == 'Circle'){
    
    for(j in 1:nrow(BinPercMat)){
      if(!BinPercMat[j,2] & !BinPercMat[j,ncol(BinPercMat)] & any(BinPercMat[j,-1])){
        BinPercMat[j,1] <- FALSE
      }
    }
    
  }
  
  
  
  for(i in 2:(ncol(BinPercMat)-1)){
    for(j in 1:nrow(BinPercMat)){
      if(!BinPercMat[j,i-1] & !BinPercMat[j,i+1] & any(BinPercMat[j,-i])){
        BinPercMat[j,i] <- FALSE
      }
    }
  }
  
  
  
  
  if(Structure == 'Circle'){
    
    for(j in 1:nrow(BinPercMat)){
      if(!BinPercMat[j,1] & !BinPercMat[j,ncol(BinPercMat)-1] & any(BinPercMat[j,-ncol(BinPercMat)])){
        BinPercMat[j,ncol(BinPercMat)] <- FALSE
      }
    }
    
  }
  
  
  
  
  
  
  if(ComputeOverlaps){
    
    OvefLapCat <- list()
    
    for(i in 1:nrow(BinPercMat)){
      
      if(i == nrow(BinPercMat)){
        if(Structure == 'Circle'){
          OvefLapCat[[i]] <- BinPercMat[i,] & BinPercMat[1,]
        }
      } else {
        OvefLapCat[[i]] <- BinPercMat[i,] & BinPercMat[i+1,]
      }
      
    }
    
    BinPercMatExt <- NULL
    
    for(i in 1:nrow(BinPercMat)){
      
      if(i == 1){
        if(Structure == 'Circle'){
          BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[1]] & !OvefLapCat[[nrow(BinPercMat)]],
                                 OvefLapCat[[i]])
        } else {
          BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[1]],
                                 OvefLapCat[[i]])
        }
        next()
      }
      
      if(i == nrow(BinPercMat)){
        if(Structure == 'Circle'){
          BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i-1]] & !OvefLapCat[[1]],
                                 OvefLapCat[[i]])
        } else {
          BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i-1]])
        }
        next()
      }
      
      BinPercMatExt <- rbind(BinPercMatExt, BinPercMat[i,] & !OvefLapCat[[i]] & !OvefLapCat[[i-1]],
                               OvefLapCat[[i]])
      
    }
    
    OverlapCat <- paste(levels(Categories)[-length(levels(Categories))],
                        levels(Categories)[-1], sep = "+")
    
    if(Structure == 'Circle'){
      OverlapCat <- c(OverlapCat,
                      paste(levels(Categories)[length(levels(Categories))],
                            levels(Categories)[1], sep = "+")
      )
    }
    
    if(Structure == 'Circle'){
      rownames(BinPercMatExt) <- as.vector(rbind(rownames(BinPercMat), OverlapCat))
    } else {
      rownames(BinPercMatExt) <- as.vector(rbind(rownames(BinPercMat), c(OverlapCat, NA)))[-nrow(BinPercMat)*2]
    }

  } else {
    BinPercMatExt <- BinPercMat
  }
  
  AllCat <- rownames(BinPercMatExt)
  
  if(Structure == "Circle"){
    ExtPath <- ExtPath[-length(ExtPath)]
    BinPercMatExt <- BinPercMatExt[, -ncol(BinPercMatExt)]
  }
  
  LowStages <- 1
  if(ComputeOverlaps){
    LowStages <- 1:2
  }
  
  while(any(BinPercMatExt[LowStages,ncol(BinPercMatExt)]) &
        all(!BinPercMatExt[setdiff(1:nrow(BinPercMatExt), LowStages),ncol(BinPercMatExt)])){
    BinPercMatExt <- cbind(BinPercMatExt[, ncol(BinPercMatExt)], BinPercMatExt[, -ncol(BinPercMatExt)])
    ExtPath <- c(ExtPath[length(ExtPath)], ExtPath[-length(ExtPath)])
  }
  
  HighStages <- nrow(BinPercMatExt)
  if(ComputeOverlaps & Structure == 'Circle'){
    HighStages <- c(nrow(BinPercMatExt)-1, nrow(BinPercMatExt))
  }
  
  while(any(BinPercMatExt[HighStages,1]) &
        all(!BinPercMatExt[setdiff(1:nrow(BinPercMatExt), HighStages),1])){
    BinPercMatExt <- cbind(BinPercMatExt[, -1], BinPercMatExt[, 1])
    ExtPath <- c(ExtPath[-1], ExtPath[1])
  }
  
  if(Structure == 'Circle'){
    ExtPath <- c(ExtPath, ExtPath[1])
    BinPercMatExt <- cbind(BinPercMatExt,
                           rep(FALSE, nrow(BinPercMatExt)))
  }
  
  Idxs <- apply(BinPercMatExt, 1, which)
  
  Bond <- lapply(Idxs[lapply(Idxs, length) > 1], range)
  
  print(1*BinPercMatExt)
  
  print(unlist(Bond))
  
  
  if(!S4Vectors::isSorted(unlist(Bond))){
    warning("Stages are not sequential!")
  }
  
  
  NodeOnGenes <- t(Nodes %*% t(RotatioMatrix))
  
  OrderedPoints <- OrderOnPath(PrinGraph = PrinGraph, Path = as.numeric(ExtPath), PointProjections = PointProjections)
  
  for(Idx in order(apply(NodeOnGenes, 1, var), decreasing = TRUE)[1:nGenes]){
    
    p <- ggplot2::ggplot(data = data.frame(x=cumsum(OrderedPoints$PathLen),
                                  y=NodeOnGenes[Idx,as.numeric(ExtPath)]),
                mapping = ggplot2::aes(x = x, y = y, color="PC")) +
      ggplot2::labs(x = "Pseudotime", y="Gene expression", title = rownames(NodeOnGenes)[Idx]) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
    
    RecCoord <- NULL
    
    for(i in 1:length(Idxs)){
      LowIds <- Idxs[[i]]-1
      LowIds[LowIds == 0] <- NA

      HiIds <- Idxs[[i]]+1
      HiIds[HiIds == length(ExtPath) + 1] <- NA
      
      LowCoord <- cumsum(OrderedPoints$PathLen)[LowIds]
      MidCoord <- cumsum(OrderedPoints$PathLen)[Idxs[[i]]]
      HighCoord <- cumsum(OrderedPoints$PathLen)[HiIds]
    
      RecCoord <- rbind(RecCoord, cbind(
        colMeans(rbind(LowCoord, MidCoord)),
        MidCoord,
        colMeans(rbind(MidCoord, HighCoord)),
        rep(names(Idxs)[i], length(MidCoord))
        )
      )

    }
    
    colnames(RecCoord) <- c("Min", "Med", "Max", "Stage")
    RecCoord <- data.frame(RecCoord)
    RecCoord$Min <- as.numeric(as.character(RecCoord$Min))
    RecCoord$Med <- as.numeric(as.character(RecCoord$Med))
    RecCoord$Max <- as.numeric(as.character(RecCoord$Max))
    
    
    RecCoord$Stage <- factor(as.character(RecCoord$Stage), levels = AllCat)
    
    RecCoord$Min[is.na(RecCoord$Min)] <- RecCoord$Med[is.na(RecCoord$Min)]
    RecCoord$Max[is.na(RecCoord$Max)] <- RecCoord$Med[is.na(RecCoord$Max)]
    
    p <- p + ggplot2::geom_rect(data = RecCoord, mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                       ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
      ggplot2::geom_point(data = data.frame(x=OrderedPoints$PositionOnPath, y=ExpData[,Idx]),
                 mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
      ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))
    
    print(p)
    
  }
  
  CellsPT <- OrderedPoints$PositionOnPath
  names(CellsPT) <-  rownames(ExpData)
  
  return(list(Structure = Structure, ExtPath = ExtPath,
              NodesPT = cumsum(OrderedPoints$PathLen),
              NodesExp = NodeOnGenes[order(apply(NodeOnGenes, 1, var), decreasing = TRUE),as.numeric(ExtPath)],
              CellsPT = CellsPT,
              CellExp = t(ExpData[,order(apply(ExpData, 1, var), decreasing = TRUE)]),
              RecCoord = RecCoord,
              StageOnNodes = Idxs
              )
         )
  
}





















################################################################################
#
# PlotsAndDistill ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param GeneExprMat 
#' @param StartSet 
#' @param Categories 
#' @param Topology 
#' @param IgnoreTail 
#' @param PlanVarLimit 
#' @param PlanVarLimitIC 
#' @param DistillThr 
#' @param Log 
#' @param StartQuant 
#' @param OutThr 
#' @param OutThrPCA 
#' @param Title 
#' @param MinProlCells 
#' @param PCACenter 
#' @param PlotDebug 
#' @param Mode 
#' @param ExtMode 
#' @param nNodes 
#' @param InitStructNodes 
#' @param DipPVThr 
#' @param PCAFilter 
#' @param PCAProjCenter 
#' @param StopCrit 
#' @param Filter 
#' @param CompareSet 
#'
#' @return
#' @export
#'
#' @examples
PlotsAndDistill <- function(GeneExprMat, StartSet, Categories, Topology = "Circle", IgnoreTail = FALSE, PlanVarLimit = .9, PlanVarLimitIC = .95,
                            DistillThr = .7, Log = TRUE, StartQuant = .5, OutThr = 3, OutThrPCA = 3, Title = '', MinProlCells = 20, PCACenter = FALSE,
                            PlotDebug = FALSE, Mode = "VarPC", ExtMode = 2, nNodes = 40, InitStructNodes = 20, DipPVThr = 1e-4, PCAFilter = TRUE,
                            PCAProjCenter = TRUE, StopCrit = .95, Filter = TRUE, QuaThr = .5, CompareSet = list(), EstProlif = "MeanPerc") {
  
  
  # Perform the initial analysis
  
  PrGraph.Initial <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = StartSet, OutThr = OutThr, nNodes = nNodes, QuaThr = QuaThr,
                                       VarThr = .99, Categories = Categories, GraphType = Topology, PlanVarLimit = PlanVarLimit,
                                       PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, InitStructNodes = InitStructNodes,
                                       MinBranDiff = 2, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                       PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = FALSE,
                                       PCAFilter = PCAFilter, OutThrPCA = OutThrPCA, PlotDebug = PlotDebug, EstProlif = EstProlif)
  
  # Distill genes
  
  DistilledGenes <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = Mode, DistillThr = DistillThr, FastReduce = FALSE,
                                    ExtMode = ExtMode, StopCrit = StopCrit, OutThr = OutThr, nNodes = nNodes, VarThr = .99, Categories = Categories,
                                    GraphType = Topology, PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE,
                                    InitStructNodes = InitStructNodes, MinBranDiff = 2, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                    PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = FALSE, PlotDebug = PlotDebug,
                                    IgnoreTail = IgnoreTail, StartQuant = StartQuant, PCAFilter = PCAFilter, OutThrPCA = OutThrPCA, EstProlif = EstProlif)
  
  # Perform the analysis on the final structure
  
  PrGraph.Final <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = DistilledGenes, OutThr = OutThr, nNodes = nNodes, 
                                     VarThr = .99, Categories = Categories, GraphType = Topology, PlanVarLimit = PlanVarLimit,
                                     PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, InitStructNodes = InitStructNodes,
                                     MinBranDiff = 2, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                     PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = FALSE, PCAFilter = PCAFilter,
                                     OutThrPCA = OutThrPCA, PlotDebug = PlotDebug, EstProlif = EstProlif)
  
  # Plot the initial projection
  ProjectOnPrincipalGraph(Nodes = PrGraph.Initial$PrinGraph$Nodes, Edges = PrGraph.Initial$PrinGraph$Edges, Points = PrGraph.Initial$Data,
                          UsedPoints = NULL, Categories = PrGraph.Initial$Categories, Title=paste(Title, "(Initial)"),
                          PCACenter = TRUE, ShowFitted = FALSE)
  
  # Plot the final projection
  ProjectOnPrincipalGraph(Nodes = PrGraph.Final$PrinGraph$Nodes, Edges = PrGraph.Final$PrinGraph$Edges, Points = PrGraph.Final$Data,
                          UsedPoints = NULL, Categories = PrGraph.Final$Categories, Title=paste(Title, "(Filtered)"),
                          PCACenter = TRUE, ShowFitted = FALSE)
  
  # Plot the difference in number
  barplot(c(length(StartSet), length(DistilledGenes)), ylab = "Number of genes", names.arg = c("Initial", "Filtered"))
  
  PTList <- list()
  
  if(length(CompareSet) > 0){
    
    CIListEnd <- list()
    CIListStart <- list()
    
    for(i in 1:length(CompareSet)){
      CIListEnd[[length(CIListEnd)+1]] <- prop.test(length(intersect(CompareSet[[i]], DistilledGenes)), length(DistilledGenes), conf.level = .95)$conf.int[1:2]*100
      CIListStart[[length(CIListStart)+1]] <- prop.test(length(intersect(CompareSet[[i]], StartSet)), length(StartSet), conf.level = .95)$conf.int[1:2]*100
    }
    
    IntersectListEnd <- lapply(CompareSet, function(x) intersect(x, DistilledGenes))
    IntersectListStart <- lapply(CompareSet, function(x) intersect(x, StartSet))
    
    for(i in 1:length(CompareSet)){
      
      PTList[[i]] <- prop.test(c(length(IntersectListEnd[[i]]), length(IntersectListStart[[i]])), c(length(DistilledGenes), length(StartSet)))
      
      yMax <- max(c(CIListEnd[[i]], CIListStart[[i]]))
      
      B <- barplot(100*c(length(IntersectListEnd[[i]]), length(IntersectListStart[[i]]))/c(length(DistilledGenes), length(StartSet)),
                   names.arg = c("Filtered", "Initial"), ylab = "Percentage of genes identified", ylim = c(0, yMax), main = names(CompareSet)[i])
      
      arrows(x0 = B[1], x1 = B[1], y0 = CIListEnd[[i]][1], y1 = CIListEnd[[i]][2], angle = 90, length = .5, lwd = 2, code = 3)
      arrows(x0 = B[2], x1 = B[2], y0 = CIListStart[[i]][1], y1 = CIListStart[[i]][2], angle = 90, length = .5, lwd = 2, code = 3)
      
    }
    
    names(PTList) <- names(CompareSet)
    
  }
  
  return(list(Genes = DistilledGenes, PTList = PTList, FinalStruct = PrGraph.Final, InitialStruct = PrGraph.Initial))
  
}


























################################################################################
#
# PlotOnStages ------------------------------------------------------
#
################################################################################


#' Title
#'
#' @param GeneExprMat 
#' @param StartSet 
#' @param Categories 
#' @param Topology 
#' @param DistillThr 
#' @param IgnoreTail 
#' @param Log 
#' @param StartQuant 
#' @param Title 
#' @param PlanVarLimit 
#' @param PlanVarLimitIC 
#' @param KeepOriginal 
#' @param PCACenter 
#' @param PlotDebug 
#' @param Mode 
#' @param ExtMode 
#' @param MaxRounds 
#' @param StopCrit 
#' @param ExpQuant 
#' @param StopMode 
#' @param Parallel 
#' @param nCores 
#' @param ClusType 
#'
#' @return
#' @export
#'
#' @examples
ProjectAndExpand <- function(GeneExprMat, StartSet, Categories, Topology = 'Circle', DistillThr = .6,
                             IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al",
                             PlanVarLimit = .85, PlanVarLimitIC = .9, KeepOriginal = TRUE, InitStructNodes = 20,
                             PCACenter = FALSE, PlotDebug = FALSE, Mode = "VarPC", ExtMode = 3,
                             MaxRounds = 15, StopCrit = .95, ExpQuant = .01, StopMode =1, EstProlif = "MeanPerc",
                             Parallel = TRUE, nCores = NULL, ClusType = "PSOCK", QuaThr = .5, MinProlCells = 20) {
  
  
  # Produce the initial analysis
  
  Info.Exp <- PlotsAndDistill(GeneExprMat = GeneExprMat, StartSet = StartSet, QuaThr = QuaThr, MinProlCells = MinProlCells,
                              Categories = Categories, Topology = Topology, DistillThr = DistillThr, InitStructNodes = InitStructNodes,
                              IgnoreTail = IgnoreTail, Log = Log, StartQuant = StartQuant, Title = Title, PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC,
                              PCACenter = PCACenter, PlotDebug = PlotDebug, Mode = Mode, ExtMode = ExtMode, EstProlif = EstProlif)
  
  # Produce the initial processed data
  
  Proc.Exp <- PlotOnStages(Structure = Topology, TaxonList = Info.Exp$FinalStruct$TaxonList[[length(Info.Exp$FinalStruct$TaxonList)]],
                           Categories = Info.Exp$FinalStruct$Categories, nGenes = 2, 
                           PrinGraph = Info.Exp$FinalStruct$PrinGraph,
                           Net = Info.Exp$FinalStruct$Net[[length(Info.Exp$FinalStruct$Net)]],
                           SelThr = .35, ComputeOverlaps = TRUE, ExpData = Info.Exp$FinalStruct$FiltExp,
                           RotatioMatrix = Info.Exp$FinalStruct$PCAData$rotation[,1:Info.Exp$FinalStruct$nDims],
                           PointProjections = Info.Exp$FinalStruct$ProjPoints[[length(Info.Exp$FinalStruct$ProjPoints)]])
  
  
  TaxonList <- Info.Exp$FinalStruct$TaxonList[[length(Info.Exp$FinalStruct$TaxonList)]]
  TaxVect <- rep(NA, ncol(Proc.Exp$NodesExp)-1)
  
  for(i in 1:length(TaxonList)){
    TaxVect[TaxonList[[i]]] <- i 
  }
  
  AllMean <- apply(log10(GeneExprMat[, colnames(Proc.Exp$CellExp)] + 1), 1, mean)
  SelMean <- apply(log10(GeneExprMat[rownames(Proc.Exp$CellExp), colnames(Proc.Exp$CellExp)] + 1), 1, mean)
  
  print(paste(sum(AllMean < min(SelMean)), "genes will be excluded a priori from the analysis"))
  print("Computing median of IQR/median")
  
  if(Parallel){
    
    tictoc::tic()
    if(is.null(nCores)){
      nCores <- parallel::detectCores() - 1
    }
    
    cl <- parallel::makeCluster(nCores, type = ClusType)
    
    parallel::clusterExport(cl, "TaxVect", environment())
    
    MedianVar <- parallel::parApply(cl,log10(GeneExprMat[AllMean >= min(SelMean), colnames(Proc.Exp$CellExp)] + 1), 1, function(x){
      AGG <- aggregate(x, by=list(TaxVect), median)
      AGGVect <- AGG[,2]
      names(AGGVect) <- paste(AGG[,1])
      median((aggregate(x, by=list(TaxVect), IQR))[,2]/AGGVect, na.rm=TRUE)
    })
    
    parallel::stopCluster(cl)
    
    tictoc::toc()
    
  } else {
    
    tictoc::tic()
    MedianVar <- apply(log10(GeneExprMat[AllMean >= min(SelMean), colnames(Proc.Exp$CellExp)] + 1), 1, function(x){
      AGG <- aggregate(x, by=list(TaxVect), median)
      AGGVect <- AGG[,2]
      names(AGGVect) <- paste(AGG[,1])
      median((aggregate(x, by=list(TaxVect), IQR))[,2]/AGGVect, na.rm=TRUE)
    })
    tictoc::toc()
    
  }
  
  
  
  MedianVar <- MedianVar[!is.infinite(MedianVar)]
  MedianVar <- MedianVar[!is.na(MedianVar)]
  
  PrevDistilled <- names(MedianVar) %in% rownames(Proc.Exp$NodesExp)
  
  boxplot(MedianVar ~ PrevDistilled)
  
  # wilcox.test(MedianVar ~ PrevDistilled)
  # t.test(MedianVar ~ PrevDistilled)
  
  GeneStageList <- list()
  Info.Exp.List <- list()
  Proc.Exp.List <- list()
  
  GeneStageList[[1]] <- unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], ExpQuant))),
                                 rownames(Proc.Exp$NodesExp)))
  
  length(GeneStageList[[1]])/nrow(Proc.Exp$NodesExp)
  
  DONE <- FALSE
  Info.Exp.Final <- NULL
  Round.Count <- 0
  
  while(!DONE){
    
    # Repeating the analysis untill convergence (or max iterations)
    
    Info.Exp.StageI <- PlotsAndDistill(GeneExprMat = GeneExprMat, StartSet = GeneStageList[[length(GeneStageList)]],
                                       Categories = Categories, Topology = Topology, DistillThr = DistillThr,
                                       IgnoreTail = IgnoreTail, Log = Log, StartQuant = StartQuant, Title = paste(Title, "STEP", Round.Count), PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC,
                                       PCACenter = PCACenter, PlotDebug = PlotDebug, Mode = Mode, ExtMode = ExtMode)
    
    Info.Exp.List[[length(Info.Exp.List)+1]] <- Info.Exp.StageI
    
    
    Proc.Exp.StageI <- PlotOnStages(Structure = Topology, TaxonList = Info.Exp.StageI$FinalStruct$TaxonList[[length(Info.Exp.StageI$FinalStruct$TaxonList)]],
                                    Categories = Info.Exp.StageI$FinalStruct$Categories, nGenes = 2,
                                    PrinGraph = Info.Exp.StageI$FinalStruct$PrinGraph,
                                    Net = Info.Exp.StageI$FinalStruct$Net[[length(Info.Exp.StageI$FinalStruct$Net)]],
                                    SelThr = .35, ComputeOverlaps = TRUE, ExpData = Info.Exp.StageI$FinalStruct$FiltExp,
                                    RotatioMatrix = Info.Exp.StageI$FinalStruct$PCAData$rotation[,1:Info.Exp.StageI$FinalStruct$nDims],
                                    PointProjections = Info.Exp.StageI$FinalStruct$ProjPoints[[length(Info.Exp.StageI$FinalStruct$ProjPoints)]])
    
    Proc.Exp.List[[length(Proc.Exp.List)+1]] <- Proc.Exp.StageI
    
    print(paste(length(Proc.Exp.StageI$NodesExp), "genes selected"))
    
    OnlyOld <- setdiff(GeneStageList[[length(GeneStageList)]], rownames(Proc.Exp.StageI$NodesExp))    
    print(paste(length(OnlyOld), "genes removed:"))
    print(OnlyOld)
    
    OnlyNew <- setdiff(rownames(Proc.Exp.StageI$NodesExp), GeneStageList[[length(GeneStageList)]])
    print(paste(length(OnlyNew), "genes added:"))
    print(OnlyNew)
    
    if(StopMode == 1){
      if( (length(OnlyOld) + length(OnlyNew))/length(GeneStageList[[length(GeneStageList)]]) < 1 - StopCrit ){
        print("Converged! - Mode 1")
        DONE <- TRUE
        break()
      }
    }
    
    Round.Count <- Round.Count + 1
    
    print(paste("Round", Round.Count, "Done"))
    
    if(Round.Count >= MaxRounds){
      print("Maximum number of iteration reached")
      DONE <- TRUE
      break()
    }
    
    TaxonList <- Info.Exp.StageI$FinalStruct$TaxonList[[length(Info.Exp.StageI$FinalStruct$TaxonList)]]
    TaxVect <- rep(NA, ncol(Proc.Exp.StageI$NodesExp)-1)
    
    for(i in 1:length(TaxonList)){
      TaxVect[TaxonList[[i]]] <- i 
    }
    
    AllMean <- apply(log10(GeneExprMat[, colnames(Proc.Exp.StageI$CellExp)] + 1), 1, mean)
    SelMean <- apply(log10(GeneExprMat[rownames(Proc.Exp.StageI$CellExp), colnames(Proc.Exp.StageI$CellExp)] + 1), 1, mean)
    
    print(paste(sum(AllMean < min(SelMean)), "genes will be excluded a priori from the analysis"))
    print("Computing median of IQR/median")
    
    
    
    
    
    
    
    if(Parallel){
      
      tictoc::tic()
      if(is.null(nCores)){
        nCores <- parallel::detectCores() - 1
      }
      
      cl <- parallel::makeCluster(nCores, type = ClusType)
      
      parallel::clusterExport(cl, "TaxVect", environment())
      
      MedianVar <- parallel::parApply(cl,log10(GeneExprMat[AllMean >= min(SelMean), colnames(Proc.Exp.StageI$CellExp)] + 1), 1, function(x){
        AGG <- aggregate(x, by=list(TaxVect), median)
        AGGVect <- AGG[,2]
        names(AGGVect) <- paste(AGG[,1])
        median((aggregate(x, by=list(TaxVect), IQR))[,2]/AGGVect, na.rm=TRUE)
      })
      
      parallel::stopCluster(cl)
      tictoc::toc()
      
    } else {
      
      tictoc::tic()
      MedianVar <- apply(log10(GeneExprMat[AllMean >= min(SelMean), colnames(Proc.Exp.StageI$CellExp)] + 1), 1, function(x){
        AGG <- aggregate(x, by=list(TaxVect), median)
        AGGVect <- AGG[,2]
        names(AGGVect) <- paste(AGG[,1])
        median((aggregate(x, by=list(TaxVect), IQR))[,2]/AGGVect, na.rm=TRUE)
      })
      tictoc::toc()
      
    }
    
    MedianVar <- MedianVar[!is.infinite(MedianVar)]
    MedianVar <- MedianVar[!is.na(MedianVar)]
    
    PrevDistilled <- names(MedianVar) %in% rownames(Proc.Exp.StageI$NodesExp)
    
    boxplot(MedianVar ~ PrevDistilled, main = paste(Title, "STEP", Round.Count))
    
    # wilcox.test(MedianVar ~ PrevDistilled)
    # t.test(MedianVar ~ PrevDistilled)
    
    if(KeepOriginal){
      GeneStageList[[length(GeneStageList)+1]] <- unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], ExpQuant))),
                                                           rownames(Proc.Exp$NodesExp)))
      GeneNumbVarRat <- length(GeneStageList[[length(GeneStageList)]])/nrow(Proc.Exp$NodesExp)
    } else {
      GeneStageList[[length(GeneStageList)+1]] <- unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], ExpQuant))),
                                                           rownames(Proc.Exp.StageI$NodesExp)))
      GeneNumbVarRat <- length(GeneStageList[[length(GeneStageList)]])/nrow(Proc.Exp.StageI$NodesExp)
    }
    
    GeneIntRatPrev <- length(intersect(GeneStageList[[length(GeneStageList)]], GeneStageList[[length(GeneStageList)-1]]))/
      length(GeneStageList[[length(GeneStageList)-1]])
    GeneIntRatCurr <- length(intersect(GeneStageList[[length(GeneStageList)]], GeneStageList[[length(GeneStageList)-1]]))/
      length(GeneStageList[[length(GeneStageList)]])
    
    print(paste("Gene number variation", GeneNumbVarRat))
    print(paste("Gene intersection variation (Prev)", GeneIntRatPrev))
    print(paste("Gene intersection variation (Curr)", GeneIntRatCurr))
    
    Sys.sleep(10)
    
    if(StopMode == 2){
      if(GeneIntRatPrev > StopCrit & GeneIntRatCurr > StopCrit){
        print("Converged! - Mode 2")
        DONE <- TRUE
        break()
      }
    }
    
  }
  
  return(list(StartInfo = Info.Exp, StartProc = Proc.Exp,
              ListInfo = Info.Exp.List, ListProc = Proc.Exp.List,
              GeneList = GeneStageList))
}














