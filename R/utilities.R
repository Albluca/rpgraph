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




#' Title
#'
#' @param x 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
CircShift <- function(x, n = 1) {
  if(n == 0){
    x
  } else {
    c(tail(x, -n), head(x, n))
  } 
}




#' Title
#'
#' @param XVect 
#' @param YVect 
#'
#' @return
#' @export
#'
#' @examples
CircCor <- function(XVect, YVect, method = "pearson") {
  
  LocShift <- function(i) {
    return(CircShift(XVect, i))
  }
  
  AllShift <- sapply(X = 0:(length(XVect)-1), LocShift)
  
  LocCor <- function(XVect) {
    return(cor(XVect, YVect, method = method))
  }
  
  return(apply(AllShift, 2, LocCor))
}




SmoothFilter <- function(CateVect, Weigth, Thr) {
  
  if(is.na(CateVect[1]) | is.na(CateVect[length(CateVect)])){
    return(CateVect)
  }
  
  if(length(unique(CateVect))==1){
    return(CateVect)
  }
  
  if(CateVect[1] == CateVect[length(CateVect)]){
    if(sum(CateVect == CateVect[1], na.rm = TRUE) > length(CateVect)/2){
      RefSizes <- Weigth[CateVect == CateVect[1]]
      TarSizes <- Weigth[CateVect != CateVect[1]]
      
      if((mean(TarSizes) - mean(RefSizes))/sd(RefSizes) < -Thr){
        CateVect[] <- CateVect[1]
      }
      return(CateVect)
    }
  }
  
  return(CateVect)
  
}




#' Title
#'
#' @param StageMatrix 
#' @param NodePenalty 
#'
#' @return
#' @export
#'
#' @examples
FitStagesCirc <- function(StageMatrix, NodePenalty, QuantSel = .05) {
  
  NormStageMatrix <- StageMatrix
  
  # Find which columns are associated with at least one stage
  
  ToAnalyze <- which(colSums(StageMatrix)>0)
  
  if(length(ToAnalyze)<nrow(StageMatrix)){
  
    # Not enough columns. Padding around them
      
    PaddedIdxs <- c(
      (min(ToAnalyze) - nrow(StageMatrix)):(min(ToAnalyze)-1),
      ToAnalyze,
      (max(ToAnalyze)+1):(max(ToAnalyze) + nrow(StageMatrix))
    )
    
    PaddedIdxs[PaddedIdxs > ncol(StageMatrix)] <-  PaddedIdxs[PaddedIdxs > ncol(StageMatrix)] - ncol(StageMatrix)
    PaddedIdxs[PaddedIdxs <= 0] <-  PaddedIdxs[PaddedIdxs <= 0] + ncol(StageMatrix)
    
    ToAnalyze <- sort(unique(PaddedIdxs))
    
  }
  
  # Selecting the matrix that will be analysed and saving the indices
  
  NormStageMatrix <- NormStageMatrix[,ToAnalyze]
  NormStageMatrixIdx <- ToAnalyze
  
  NormNodePenalty <- NodePenalty[ToAnalyze]
  
  Possibilities <- combn(1:ncol(NormStageMatrix), nrow(NormStageMatrix))
  
  ToKeep <- !apply(Possibilities, 2, is.unsorted)
  
  Possibilities <- Possibilities[, ToKeep]
  dim(Possibilities) <- c(length(Possibilities)/length(ToKeep), length(ToKeep))
  
  PathPenelity <- function(ChangeNodes, InitialStage) {
    
    Sphases <- rep(InitialStage, ncol(NormStageMatrix))
    
    for (i in 1:(length(ChangeNodes)-1)) {
      Sphases[ChangeNodes[i]:ChangeNodes[length(ChangeNodes)]] <- InitialStage + i
    }
    
    Sphases[Sphases>length(ChangeNodes)] <- Sphases[Sphases>length(ChangeNodes)] - length(ChangeNodes)
    
    sum(NormNodePenalty*(apply(NormStageMatrix, 2, max) - mapply("[[", apply(NormStageMatrix, 2, as.list), Sphases))^2)
    
  }
  
  CombinedInfo <- NULL
  
  for (i in 1: nrow(NormStageMatrix)) {
    Penalities <- apply(Possibilities, 2, PathPenelity, InitialStage = i)
    Best <- which(Penalities == min(Penalities))
    BestPossibilities <- Possibilities[, Best]
    dim(BestPossibilities) <- c(nrow(NormStageMatrix), length(BestPossibilities)/nrow(NormStageMatrix))
    
    CombinedInfo <- cbind(CombinedInfo, rbind(rep(i, ncol(BestPossibilities)),
                                              rep(min(Penalities), ncol(BestPossibilities)),
                                              BestPossibilities)
    )
  }
  
  SelIdxs <- which(CombinedInfo[2,] <= quantile(CombinedInfo[2,], QuantSel))
  
  SelInfo <- CombinedInfo[,SelIdxs]
  dim(SelInfo) <- c(length(SelInfo)/length(SelIdxs), length(SelIdxs))

  ChangeNodes <- NormStageMatrixIdx[SelInfo[-c(1:2),]]
  dim(ChangeNodes) <- c(nrow(SelInfo)-2, length(ChangeNodes)/(nrow(SelInfo)-2))
  
  ExpandStages <- function(idx) {
    
    StageVect <- rep(SelInfo[1,idx], ncol(StageMatrix))
    
    for (i in 1:(nrow(ChangeNodes)-1)) {
      StageVect[ChangeNodes[i,idx]:ChangeNodes[nrow(ChangeNodes),idx]] <- SelInfo[1,idx] + i
    }
    
    StageVect[StageVect>nrow(ChangeNodes)] <- StageVect[StageVect>nrow(ChangeNodes)] - nrow(ChangeNodes)
    
    return(StageVect)
    
  }
  
  return(list(Order = lapply(as.list(1:ncol(ChangeNodes)), ExpandStages), Penality = SelInfo[2,]))
  
}





#' Title
#'
#' @param StageAssociation 
#'
#' @return
#' @export
#'
#' @examples
ExtendImplicitStaging <- function(StageAssociation) {

  AllGenes <- NULL
    
  for (Stage in 1:length(StageAssociation$Stages)) {
    
    if(exists(paste("S", Stage, "_U", sep = ""), where=StageAssociation)) {
      AllGenes <- c(AllGenes, unlist(StageAssociation[paste("S", Stage, "_U", sep = "")], use.names = FALSE))
    }
    
    if(exists(paste("S", Stage, "_D", sep = ""), where=StageAssociation)) {
      AllGenes <- c(AllGenes, unlist(StageAssociation[paste("S", Stage, "_D", sep = "")], use.names = FALSE))
    }
    
  }
  
  AllGenes <- unique(AllGenes)
  
  ContributionMatrix <- NULL
  NameVect <- NULL
  
  for (Stage in 1:length(StageAssociation$Stages)) {
    
    Name <- paste("S", Stage, "_U", sep = "")
    
    if(exists(Name, where=StageAssociation)) {
      
      NameVect <- c(NameVect, Name)
      FoundGenes <- unlist(StageAssociation[Name], use.names = FALSE)
      
      ContributionMatrix <- rbind(ContributionMatrix, AllGenes %in% FoundGenes)
      
    } else {
      
      NameVect <- c(NameVect, Name)
      FoundGenes <- ''
      
      ContributionMatrix <- rbind(ContributionMatrix, AllGenes %in% FoundGenes)
    }
    
    Name <- paste("S", Stage, "_D", sep = "")
    
    if(exists(Name, where=StageAssociation)) {
      
      NameVect <- c(NameVect, Name)
      FoundGenes <- unlist(StageAssociation[Name], use.names = FALSE)
      
      ContributionMatrix <- rbind(ContributionMatrix, AllGenes %in% FoundGenes)
      
    } else {
      
      NameVect <- c(NameVect, Name)
      FoundGenes <- ''
      
      ContributionMatrix <- rbind(ContributionMatrix, AllGenes %in% FoundGenes)
    }
    
  }
  
  colnames(ContributionMatrix) <- AllGenes
  rownames(ContributionMatrix) <- NameVect
  
  UpIdx <- seq(1, nrow(ContributionMatrix), 2)
  DownIdx <- seq(2, nrow(ContributionMatrix), 2)
  
  for (i in 1:ncol(ContributionMatrix)) {
    
      if(any(ContributionMatrix[UpIdx, i])){
        ContributionMatrix[DownIdx, i] <- !ContributionMatrix[UpIdx, i]
        next()
      }
      
    if(any(ContributionMatrix[DownIdx, i])){
      ContributionMatrix[UpIdx, i] <- !ContributionMatrix[DownIdx, i]
      next()
    }
  }
  
  if(any(colSums(ContributionMatrix) != colSums(!ContributionMatrix))){
    warning("Something went orribly wrong ...")
  }
  
  for(i in 1:nrow(ContributionMatrix)){
    StageAssociation[[rownames(ContributionMatrix)[i]]] <- colnames(ContributionMatrix)[ContributionMatrix[i, ]]
  }
  
  return(StageAssociation)
 
}

