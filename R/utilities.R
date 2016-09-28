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



