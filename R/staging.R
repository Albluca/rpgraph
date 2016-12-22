StagingByGenes <- function(StageAssociation, ExpressionMatrix, NormExpressionMatrix,
                           NodeOnGenesOnPath, UsedPath, NodeSize, NodePower, LowQ, TopQ, MaxExp, MinExp,
                           nPoints, MinWit, PercNorm, StagingMode, CutOffVar) {
  
  
  print("Gene/Stage information found. Trying to Optimize")
  
  StageMatU <- NULL
  StageMatD <- NULL
  GeneCount <- rep(0, length(StageAssociation$Stages))
  
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
        
        SignStageMat <- sign(t(StageTracks) - apply(StageTracks, 2, quantile, TopQ)) * sign(t(StageTracks) - MinExp)
        
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
        
        SignStageMat <- sign(t(StageTracks) - apply(StageTracks, 2, quantile, LowQ)) * sign(t(StageTracks) - MaxExp)
        
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
  
  NoNormSummaryStageMat <- SummaryStageMat
  NoNormWeigth <- NodeSize^NodePower
  
  SummaryStageMat[SummaryStageMat<MinWit] <- 0
  
  if(PercNorm){
    SummaryStageMat <- SummaryStageMat/GeneCount
  }
  
  SummaryStageMat[is.nan(SummaryStageMat)] <- 0
  
  print("The following staging matrix will be used")
  colnames(SummaryStageMat) <- UsedPath
  print(SummaryStageMat)
  
  tictoc::tic()
  print("Direct staging")
  Staging <- FitStagesCirc(StageMatrix = SummaryStageMat,
                           NodePenalty = NodeSize^NodePower,
                           Mode = StagingMode)
  tictoc::toc()
  
  tictoc::tic()
  print("Reverse staging")
  StagingRev <- FitStagesCirc(StageMatrix = SummaryStageMat[, rev(1:ncol(SummaryStageMat))],
                              NodePenalty = rev(NodeSize)^NodePower,
                              Mode = StagingMode)  
  tictoc::toc()
  
  AllPenality <- rbind(cbind(Staging$Penality, StagingRev$Penality), rep(1:2, each=ncol(Staging$Penality)))
  
  # Idxs <- sample(x = 1:ncol(AllPenality), size = nStages, prob = max(AllPenality[2, ]) - AllPenality[2, ])
  # Idxs <- unique(Idxs)
  
  # Idxs <- order(AllPenality[2, ])[1:nStages]
  
  Idxs <- which(AllPenality[2, ] == min(AllPenality[2, ]))
  
  SelPenality <- AllPenality[,Idxs]
  dim(SelPenality) <- c(4, length(SelPenality)/4)
  
  DirectPenality <- NULL
  DirectChanges <- NULL
  
  if(sum(SelPenality[4,] == 1) > 0){
    
    
    ExpandStages <- function(idx) {
      
      ChangeNodes <- Staging$Possibilities[ , SelPenality[3, idx]]
      
      StageVect <- rep(SelPenality[1, idx], ncol(SummaryStageMat))
      
      tStart <- NA
      tEnd <- NA
      
      for (i in 1:(length(ChangeNodes)-1)) {
        
        if(!is.na(ChangeNodes[i])){
          tStart <- ChangeNodes[i]
        }
        
        tEnd <- ChangeNodes[i+1]
        
        if(is.na(tStart) | is.na(tEnd)){
          next()
        }
        
        StageVect[tStart:(tEnd-1)] <- SelPenality[1, idx] + i
      }
      
      StageVect[StageVect>length(ChangeNodes)] <- StageVect[StageVect>length(ChangeNodes)] - length(ChangeNodes)
      
      return(StageVect)
      
    }
    
    SelPenIdx <- which(SelPenality[4,]==1)
    
    DirectPenality <- SelPenality[2,SelPenIdx]
    DirectChanges <- t(sapply(SelPenIdx, ExpandStages))
    
  }
  
  ReversePenality <- NULL
  ReverseChanges <- NULL
  
  if(sum(SelPenality[4,] == 2) > 0){
    
    ExpandStages <- function(idx) {
      
      ChangeNodes <- StagingRev$Possibilities[ , SelPenality[3, idx]]
      
      StageVect <- rep(SelPenality[1, idx], ncol(SummaryStageMat))
      
      tStart <- NA
      tEnd <- NA
      
      for (i in 1:(length(ChangeNodes)-1)) {
        
        if(!is.na(ChangeNodes[i])){
          tStart <- ChangeNodes[i]
        }
        
        tEnd <- ChangeNodes[i+1]
        
        if(is.na(tStart) | is.na(tEnd)){
          next()
        }
        
        StageVect[tStart:(tEnd-1)] <- SelPenality[1, idx] + i
      }
      
      StageVect[StageVect>length(ChangeNodes)] <- StageVect[StageVect>length(ChangeNodes)] - length(ChangeNodes)
      
      return(rev(StageVect))
      
    }
    
    SelPenIdx <- which(SelPenality[4,]==2)
    
    ReversePenality <- SelPenality[2,SelPenIdx]
    ReverseChanges <- t(sapply(SelPenIdx, ExpandStages))
    
  }
  
  AllStg <- rbind(DirectChanges, ReverseChanges)
  AllPen <- c(DirectPenality, ReversePenality)
  AllDir <- c(rep("Dir", length(DirectPenality)),
              rep("Rev", length(ReversePenality)))
  colnames(AllStg) <- UsedPath
  
  return(list(AllStg = AllStg, AllPen = AllPen, AllDir = AllDir,
              SummaryStageMat = SummaryStageMat, NoNormWeigth = NoNormWeigth))
  
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
FitStagesCirc <- function(StageMatrix, NodePenalty, Mode = 1) {
  
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
  
  Possibilities <- combn(c(1:ncol(NormStageMatrix), rep(NA, nrow(NormStageMatrix))), nrow(NormStageMatrix))
  
  NoChange <- apply(is.na(Possibilities), 2, sum) == nrow(NormStageMatrix)
  
  ToKeep <- (!apply(Possibilities, 2, is.unsorted, na.rm = TRUE)) &
    (apply(!is.na(Possibilities), 2, sum) != 1) &
    (!is.na(Possibilities[nrow(NormStageMatrix),]))
    
  Possibilities <- Possibilities[, ToKeep | NoChange]
  dim(Possibilities) <- c(length(Possibilities)/(sum(ToKeep)+1), sum(ToKeep)+1)
  
  PathPenality <- function(ChangeNodes, InitialStage, Mode) {
    
    Sphases <- rep(InitialStage, ncol(NormStageMatrix))
    
    tStart <- NA
    tEnd <- NA
    
    for (i in 1:(length(ChangeNodes)-1)) {
      
      if(!is.na(ChangeNodes[i])){
        tStart <- ChangeNodes[i]
      }
      
      tEnd <- ChangeNodes[i+1]
      
      if(is.na(tStart) | is.na(tEnd)){
        next()
      }
      
      Sphases[tStart:(tEnd-1)] <- InitialStage + i
    }
    
    Sphases[Sphases>length(ChangeNodes)] <- Sphases[Sphases>length(ChangeNodes)] - length(ChangeNodes)
    
    if(Mode == 1){
      # Squared distance from the maximum
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, sum) - diag(NormStageMatrix[Sphases,]))
        )
      )
    }
    
    if(Mode == 2){
      # Sum of "off-stage" contributions
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, max) - mapply("[[", apply(NormStageMatrix, 2, as.list), Sphases))
        )
      )
    }
    
    if(Mode == 3){
      # Sum binary difference from maximum
      return(
        sum(
          NormNodePenalty*(apply(NormStageMatrix, 2, which.max) != Sphases)
        )
      )
    }
    
    if(Mode == 4){
      # Sum binary difference from maximum
      tDiff <- abs(apply(NormStageMatrix, 2, which.max) - Sphases)
      tDiff[tDiff > nrow(NormStageMatrix)/2] <- nrow(NormStageMatrix) - tDiff[tDiff > nrow(NormStageMatrix)/2]
      return(
        sum(
          NormNodePenalty*tDiff
        )
      )
    }
    
  }
  
  CombinedInfo <- NULL
  
  for (i in 1:nrow(NormStageMatrix)) {
    CombinedInfo <- cbind(CombinedInfo,
                          rbind(rep(i, sum(ToKeep)+1),
                                apply(Possibilities, 2, PathPenality, InitialStage = i, Mode = Mode),
                                1:(sum(ToKeep)+1)
                          )
    )
  }
  
  return(list(Penality = CombinedInfo, Possibilities = Possibilities))
  
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


