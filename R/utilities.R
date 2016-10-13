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




#' Title
#'
#' @param SETVM 
#'
#' @return
#'
#' @examples
SelectVM <- function(SETVM = NULL) {
  
  # Check operating Systems and set variable if necessary
  
  BaseJavaDir <- "/Library/Java/JavaVirtualMachines/"
  
  if(grep("apple", R.version$platform)){
    # R has been compiled for Mac
    print("Mac OS detected. Setting up Java environmental variables.")
    
    # Get environmental variables
    
    if(Sys.getenv("JAVA_HOME") == ""){
      
      AvailableVMS <- list.dirs(BaseJavaDir, recursive = FALSE, full.names = FALSE)
      AvailableVMS <- AvailableVMS[grep("jdk", AvailableVMS)]
      
      print(paste(length(AvailableVMS), "JVMs found"))
      
      JVMVersionNumb <- NULL
      for(i in 1:length(AvailableVMS)){
        JVMVersionNumb <- c(JVMVersionNumb,
                            regmatches(AvailableVMS[i],
                                       regexpr("[0-9]+[.]+[0-9]+[.]+[0-9]+[_]*[0-9]*", AvailableVMS[i])
                            )
        )
      }
      
      JVM_VersionNumbType <- as.numeric_version(sub("_", ".", JVMVersionNumb, fixed = TRUE))
      
      if(!is.null(SETVM)){
        
        SelectedVMs <- AvailableVMS[grep(SETVM, AvailableVMS)]
        
        if(length(SelectedVMs)==1){
          print(paste("Using", SelectedVMs))
          
          options("java.home"=paste(BaseJavaDir, SelectedVMs, "/Contents/Home/jre", sep=""))
          dyn.load(paste(BaseJavaDir, SelectedVMs, "/Contents/Home/jre/lib/server/libjvm.dylib", sep=""))
        } else {
          VersionToUse <- max(JVM_VersionNumbType[grep(SETVM, AvailableVMS)])
          DirToUse <- AvailableVMS[which(JVM_VersionNumbType == VersionToUse)]
          
          print(paste("Using", DirToUse))
          
          options("java.home"=paste(BaseJavaDir, DirToUse, "/Contents/Home/jre", sep=""))
          dyn.load(paste(BaseJavaDir, DirToUse, "/Contents/Home/jre/lib/server/libjvm.dylib", sep=""))
        }
        
      } else {
        
        VersionToUse <- max(JVM_VersionNumbType)
        DirToUse <- AvailableVMS[which(JVM_VersionNumbType == VersionToUse)]
        
        print(paste("Using", DirToUse))
        
        options("java.home"=paste(BaseJavaDir, DirToUse, "/Contents/Home/jre", sep=""))
        dyn.load(paste(BaseJavaDir, DirToUse, "/Contents/Home/jre/lib/server/libjvm.dylib", sep=""))
      }
    }
  }
  
}



#' Sample correlation with partial order
#'
#' @param OrdVect 
#' @param Vect1 
#' @param Vect2 
#' @param Round 
#'
#' @return
#' @export
#'
#' @examples
Rand.POrd.Cor <- function(OrdVect, Vect1, Vect2, Round) {
  
  OrdVect <- factor(OrdVect)
    
  Vals <- unique(OrdVect)
  Multiplicity <- table(OrdVect)
  
  if(max(Multiplicity)==1){
    return(cor(Vect1, Vect2))
  }
  
  CorVect <- NULL
  
  for(i in 1:Round){
    
    SplitV1 <- split(Vect1, OrdVect)
    SplitV2 <- split(Vect2, OrdVect)
    
    Sampled1 <- lapply(SplitV1[lapply(SplitV1, length) > 1], sample)
    SplitV1[names(Sampled1)] <- Sampled1
    
    Sampled2 <- lapply(SplitV2[lapply(SplitV2, length) > 1], sample)
    SplitV2[names(Sampled2)] <- Sampled1
    
    CorVect <- c(CorVect, cor(unlist(SampledV1), unlist(SampledV2)))
    
  }
  
  return(CorVect)
  
}






