library(dplyr)
library(ggplot2)
library(rpgraph)
library(readr)
library(readxl)
library(Pricecycle)

BaseDir <- "~/Google Drive/Datasets/Trapnell et al - Human skeletal muscle myoblasts/"

Myoblast_d0 <- read.delim(paste(BaseDir, "Myoblast_d0.txt", sep=''), stringsAsFactors=FALSE, row.names = 1)
Myoblast_d0 <- Myoblast_d0[, colnames(Myoblast_d0) != "X"]

Myoblast <- read.delim(paste(BaseDir, "Myoblast.txt", sep=''), stringsAsFactors=FALSE, row.names = 1)
Myoblast <- Myoblast[, colnames(Myoblast) != "X"]

Myoblast_ID <-  read.delim(paste(BaseDir, "Trapnell_sample_ID_SC.txt", sep=''), stringsAsFactors=FALSE)

AllMyoBlast_Raw <- t(Myoblast)
GroupName <- rownames(AllMyoBlast_Raw)

table(Myoblast_ID$GroupID)

# Data will not be log trasnformed

G1S_genes <- unlist(union(StageAssociation_Whit_G0$S2_U, StageAssociation_Whit_G0$S3_U))
G2M_genes <- unlist(union(StageAssociation_Whit_G0$S4_U, StageAssociation_Whit_G0$S5_U))

SelectedSamples <- AllMyoBlast_Raw

SampVect <- Myoblast_ID$GroupID
names(SampVect) <- Myoblast_ID$SampleID

SampVect[which(SampVect == "Myoblast_d0")]

Categories <- SampVect[rownames(SelectedSamples[names(which(SampVect == "Myoblast_d0")), ])]
AllData <- t(data.matrix(SelectedSamples[names(which(SampVect == "Myoblast_d0")), ]))

# rownames(AllData) <- toupper(Murine_HSC_Data$lp)
# AllData <- t(AllData)
# AllData[is.na(AllData)] <- 0


Output.Trap_D0 <- SelectGenesOnGraph(DataSet = AllData,
                                  StartSet = HumanGenes_GOCellCycle,
                                  Categories = factor(Categories),
                                  PCACenter = TRUE, PCAProjCenter = TRUE)


sum(toupper(G0.Sel) %in% Output.Trap_D0$Genes[[length(Output.Trap_D0$Genes)]])

sum(toupper(G1.Sel) %in% Output.Trap_D0$Genes[[length(Output.Trap_D0$Genes)]])

sum(toupper(S.Sel) %in% Output.Trap_D0$Genes[[length(Output.Trap_D0$Genes)]])

sum(toupper(G2.Sel) %in% Output.Trap_D0$Genes[[length(Output.Trap_D0$Genes)]])








TargetStruct <- Output.Trap_D0$PGStructs[[length(Output.Trap_D0$PGStructs)]]
Proc.Exp.Trap <- PlotOnStages(Structure = "Circle",
                              Categories = TargetStruct$Categories,
                              nGenes = 2,
                              TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                              PrinGraph = TargetStruct$PrinGraph,
                              Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                              SelThr = .3,
                              ComputeOverlaps = TRUE,
                              ExpData = TargetStruct$FiltExp,
                              RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                              PCACenter = TargetStruct$PCAData$center,
                              PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                              OrderOnCat = FALSE,
                              SmoothPoints = 2, MinCellPerNode = 2)



pheatmap::pheatmap(Proc.Exp.Trap$NodesExp[rownames(Proc.Exp.Trap$NodesExp) %in% toupper(G1.Sel) |
                                            rownames(Proc.Exp.Trap$NodesExp) %in% toupper(S.Sel) |
                                            rownames(Proc.Exp.Trap$NodesExp) %in% toupper(G2.Sel)],
                   cluster_cols = FALSE)


plot(Proc.Exp.Trap$NodesExp["CDK4",])




Output.Trap_D0$PGStructs[[5]]$PrinGraph


ProcStruct <- Proc.Exp.Trap


GetGenesWithPeaks <- function(DataStruct, ProcStruct, FiltMax = .1, Thr = .9) {

  NormExp <- ProcStruct$NodesExp
  NormExp[NormExp < FiltMax] <- 0

  NormExp <- NormExp[rowSums(NormExp==0) < ncol(NormExp), ]

  pheatmap::pheatmap(NormExp.Bin*1, cluster_cols = FALSE)

  NormExp.Bin <- NormExp > apply(NormExp, 1, min) + (apply(NormExp, 1, max) - apply(NormExp, 1, min))*Thr

  pheatmap::pheatmap(1*NormExp.Bin[order(rowMeans(t(t(NormExp.Bin)*1:ncol(NormExp.Bin)))), ],
                     cluster_cols = FALSE, cluster_rows = FALSE)

  SinglePk_Up <- lapply(apply(NormExp.Bin, 1, which), function(x){
    if(all(min(x):max(x) %in% x)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  SinglePk_Down <- lapply(apply(!NormExp.Bin, 1, which), function(x){
    if(all(min(x):max(x) %in% x)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  Selected <- c(names(unlist(SinglePk_Up)[unlist(SinglePk_Up)]),
                names(unlist(SinglePk_Down)[unlist(SinglePk_Down)]))

  # PlotMat <- 1*NormExp.Bin[Selected,]

  PlotMat <- 1*NormExp.Bin

  PlotMat.Numb <- t(t(PlotMat)*1:ncol(PlotMat))
  PlotMat.Numb[PlotMat.Numb == 0] <- NA

  pheatmap::pheatmap(PlotMat[order(apply(PlotMat.Numb, 1, median, na.rm=TRUE)),], show_colnames = FALSE,
                     show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)

  barplot(colSums(NormExp.Bin[Selected,]))

  InfoData <- rbind(
    colSums(PlotMat[rownames(PlotMat) %in% toupper(G1.Sel),])/sum(rownames(PlotMat) %in% toupper(G1.Sel)),
    colSums(PlotMat[rownames(PlotMat) %in% toupper(S.Sel),])/sum(rownames(PlotMat) %in% toupper(S.Sel)),
    colSums(PlotMat[rownames(PlotMat) %in% toupper(G2.Sel),]/sum(rownames(PlotMat) %in% toupper(G2.Sel)))
  )


  1*(InfoData > .2)

  apply(InfoData, 2, which.max)

  barplot(InfoData[-1,], beside = TRUE)








  StageGenes.Names <- lapply(StageGenes, function(x){
    names(which(x))
  })

  barplot(unlist(lapply(StageGenes.Names, length)), las = 2, horiz = FALSE)

  return(StageGenes.Names)

}


















































#
# DataSet = AllData
# StartSet = HumanGenes_GOCellCycle
# Categories = factor(Categories)
#
# VarThr = .99
# nNodes = 40
# Log = TRUE
# Filter = TRUE
# OutThr = 5
# PCAFilter = TRUE
# OutThrPCA = 3
# GraphType = "Circle"
# PlanVarLimit = .85
# PlanVarLimitIC = .9
# MinBranDiff = 3
# InitStructNodes = 20
# ForceLasso = FALSE
# EstProlif = "MeanPerc"
# QuaThr = .5
# NonG0Cell = NULL
# DipPVThr = 1e-4
# MinProlCells = 50
# PCACenter = TRUE
# PCAProjCenter = TRUE
# PlotDebug = FALSE
# PlotIntermediate = FALSE
# AddGenePerc = 5
# SelThr1 = .95
# SelThr2 = .99
# MadsThr =  1




