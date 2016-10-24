# A function to standardize the Initial steps of the analysis

StudyCellCycles <- function(ExpressionMatrix, Grouping, GeneSet = NULL,
                            StageAssociation = NULL,
                            PathOpt = "switch", GeneOpt = 10,
                            GeneDetectedFilter = 2.5, GeneCountFilter = 2.5,
                            MinCellExp = 1, VarFilter = 0, LogTranform = TRUE, Centering = FALSE,
                            Scaling = FALSE, nDim = NULL, nPoints = 20) {
  
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
  
  print("Plotting variance distribution (For reference only ATM)")
    
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
  
  # Curve Fit ---------------------------------------------------------------
  
  print("Stage V - Path Selection")
  
  # Start be selecting the first available path.
  
  Pattern <- igraph::graph.ring(n = nPoints, directed = FALSE, mutual = FALSE, circular = FALSE)
  
  PossiblePaths <- igraph::graph.get.subisomorphisms.vf2(graph1 = Net, graph2 = Pattern)
  
  UsedPath <- PossiblePaths[[1]]$name
  
  print("Using")
  print(UsedPath)
  
  NodeOnGenes <- t(t(Results$Nodes %*% t(NormExpressionMatrixPCA$Comp)))
  
  NumericPath <- as.numeric(unlist(lapply(strsplit(UsedPath, "V_"), "[[", 2)))
  
  # if(Circular){
  #   NumericPath <- c(NumericPath, NumericPath[1])
  # }
  
  # PathProjection <- OrderOnPath(PrinGraph = Results, Path = NumericPath, PointProjections = ProjPoints)

  NodeOnGenesOnPath <- NodeOnGenes[NumericPath,]
  
  # Select path that optimize a given behaviour
  # GraphVarSorted <- sort(apply(NodeOnGenes, 2, var), index.return=TRUE)
  # sign(t(NodeOnGenesOnPath) - apply(NodeOnGenesOnPath, 2, mean))["CDC6",]
  
  
  StageAssociation <- list(Stages = c("G1", "S", "G2"),
                           S1_U = c("E2F5", "CCNE1", "CCNE2", "CDC25A", "CDC45", "CDC6",
                                    "CDKN1A", "CDKN3", "E2F1", "MCM2", "MCM6", "NPAT",
                                    "PCNA", "SLBP"),
                           S2_U = c("BRCA1", "BRCA2", "CCNG2", "CDKN2C", "DHFR",
                                    "MSH2", "NASP", "RRM1", "RRM2", "TYMS"),
                           S3_U = c("CCNA2", "CCNF", "CENPF", "TOP2A", "BIRC5", "BUB1",
                                    "BUB1B", "CCNB1", "CCNB2", "CDK1", "CDC20", "CDC25B",
                                    "CDC25C", "CDKN2D", "CENPA", "CKS1B", "CKS2", "PLK1",
                                    "AURKA", "RACGAP1", "KIF20A"))
  
  if(is.list(StageAssociation)){
    
    StageMatU <- NULL
    StageMatD <- NULL
    
    for (Stage in 1:length(StageAssociation$Stages)) {
      
      if(exists(paste("S", Stage, "_U", sep = ""), where=StageAssociation)) {
        
        StageGenes <- unlist(StageAssociation[paste("S", Stage, "_U", sep = "")], use.names = FALSE)
        
        StageTracks <- NodeOnGenesOnPath[, StageGenes]
        SignStageMat <- sign(t(StageTracks) - apply(StageTracks, 2, quantile, .90))
        
        SignStageMat[SignStageMat <= 0] <- NA
        SignStageMat[SignStageMat > 0] <- Stage
        
        StageMatU <- rbind(StageMatU, SignStageMat)
        
      }
      
      if(exists(paste("S", Stage, "_D", sep = ""), where=StageAssociation)) {
        
        StageGenes <- unlist(StageAssociation[paste("S", Stage, "_D", sep = "")], use.names = FALSE)
        
        StageTracks <- NodeOnGenesOnPath[, StageGenes]
        SignStageMat <- sign(t(StageTracks) - apply(StageTracks, 2, quantile, .90))
        
        SignStageMat[SignStageMat <= 0] <- NA
        SignStageMat[SignStageMat > 0] <- Stage
        
        StageMatD <- rbind(StageMatD, SignStageMat)
        
      }
      
    }
    
    SummaryStageMat <- NULL
    
    for (Stage in names(StageAssociation)) {
      SummaryStageMat <- rbind(SummaryStageMat, colSums(StageMat == Stage, na.rm = TRUE))
    }
    
    SummaryStageMat <- SummaryStageMat/unlist(lapply(StageAssociation, length))
    
    rownames(SummaryStageMat) <- names(StageAssociation)
    
    StageNodeAssociation <- apply(SummaryStageMat, 2, which.max)
    
    StageNodeAssociation
    
  }
  
  print("Gene/Stage information found. Trying to Optimize")
  
  # Start be selecting the first available path.
  
  
  
}


