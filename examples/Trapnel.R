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

Trap_D0.Data <- list(Analysis = Output.Trap_D0,
                   ExpMat = log10(AllData+1),
                   Cats = Categories)


#
# Factors <- scran::computeSumFactors(x = data.matrix(AllData),
#                                     sizes=c(10, 20, 30, 40),
#                                     positive = TRUE)
#
# AllData_Norm <- t(t(AllData)/Factors)
# AllData_Norm <- AllData_Norm[, Factors > 0]
#
#
# Output.Trap_D0.Norm <- SelectGenesOnGraph(DataSet = AllData_Norm,
#                                           StartSet = HumanGenes_GOCellCycle,
#                                           Categories = factor(Categories),
#                                           PCACenter = TRUE, PCAProjCenter = TRUE)
#



TargetStruct <- Trap_D0.Data$Analysis$PGStructs[[length(Trap_D0.Data$Analysis$PGStructs)]]
Proc.Trap_D0 <- PlotOnStages(Structure = "Circle",
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
                              SmoothPoints = 2,
                              MinCellPerNode = 2,
                              Title = 'Trapnell et al (D0)')

SelStageInfo.Human <- SelStageInfo
for(i in 1:length(SelStageInfo.Human)){
  SelStageInfo.Human[[i]] <- toupper(SelStageInfo.Human[[i]])
}


Staged.Trapnel.D0 <- StageWithPeaks(DataStruct = Trap_D0.Data,
                              ProcStruct = Proc.Trap_D0,
                              ComputeG0 = TRUE,
                              FiltMax = 0,
                              Thr = .8,
                              QuantSel = .8,
                              StageInfo = SelStageInfo.Human,
                              MinNodes = 2,
                              Mode = 1,
                              G0Level = 1.1,
                              Title = 'Trapnell et al D0')










Categories <- SampVect[rownames(SelectedSamples[names(which(SampVect == "Myoblast_d1")), ])]
AllData <- t(data.matrix(SelectedSamples[names(which(SampVect == "Myoblast_d1")), ]))

# rownames(AllData) <- toupper(Murine_HSC_Data$lp)
# AllData <- t(AllData)
# AllData[is.na(AllData)] <- 0


Output.Trap_D1 <- SelectGenesOnGraph(DataSet = AllData,
                                     StartSet = HumanGenes_GOCellCycle,
                                     Categories = factor(Categories),
                                     PCACenter = TRUE, PCAProjCenter = TRUE)

Trap_D1.Data <- list(Analysis = Output.Trap_D1,
                     ExpMat = log10(AllData+1),
                     Cats = Categories)



TargetStruct <- Trap_D1.Data$Analysis$PGStructs[[length(Trap_D1.Data$Analysis$PGStructs)]]
Proc.Trap_D1 <- PlotOnStages(Structure = "Circle",
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
                             SmoothPoints = 2,
                             MinCellPerNode = 2,
                             Title = 'Trapnell et al (D1)')


Staged.Trapnel.D1 <- StageWithPeaks(DataStruct = Trap_D1.Data,
                                    ProcStruct = Proc.Trap_D1,
                                    ComputeG0 = TRUE,
                                    FiltMax = 0,
                                    Thr = .8,
                                    QuantSel = .8,
                                    StageInfo = SelStageInfo.Human,
                                    MinNodes = 2,
                                    Mode = 1,
                                    G0Level = 1.1,
                                    Title = 'Trapnell et al D1')






Categories <- SampVect[rownames(SelectedSamples[names(which(SampVect == "Myoblast_d2")), ])]
AllData <- t(data.matrix(SelectedSamples[names(which(SampVect == "Myoblast_d2")), ]))


Output.Trap_D2 <- SelectGenesOnGraph(DataSet = AllData,
                                     StartSet = HumanGenes_GOCellCycle,
                                     Categories = factor(Categories),
                                     PCACenter = TRUE, PCAProjCenter = TRUE)

Trap_D2.Data <- list(Analysis = Output.Trap_D2,
                     ExpMat = log10(AllData+1),
                     Cats = Categories)



TargetStruct <- Trap_D2.Data$Analysis$PGStructs[[length(Trap_D2.Data$Analysis$PGStructs)]]
Proc.Trap_D2 <- PlotOnStages(Structure = "Circle",
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
                             SmoothPoints = 2,
                             MinCellPerNode = 2,
                             Title = 'Trapnell et al (D2)')


Staged.Trapnel.D2 <- StageWithPeaks(DataStruct = Trap_D2.Data,
                                    ProcStruct = Proc.Trap_D2,
                                    ComputeG0 = TRUE,
                                    FiltMax = 0,
                                    Thr = .8,
                                    QuantSel = .8,
                                    StageInfo = SelStageInfo.Human,
                                    MinNodes = 2,
                                    Mode = 1,
                                    G0Level = 1.1,
                                    Title = 'Trapnell et al D2')









Categories <- SampVect[rownames(SelectedSamples[names(which(SampVect == "Myoblast_d3")), ])]
AllData <- t(data.matrix(SelectedSamples[names(which(SampVect == "Myoblast_d3")), ]))



Output.Trap_D3 <- SelectGenesOnGraph(DataSet = AllData,
                                     StartSet = HumanGenes_GOCellCycle,
                                     Categories = factor(Categories),
                                     PCACenter = TRUE, PCAProjCenter = TRUE)

Trap_D3.Data <- list(Analysis = Output.Trap_D3,
                     ExpMat = log10(AllData+1),
                     Cats = Categories)



TargetStruct <- Trap_D3.Data$Analysis$PGStructs[[length(Trap_D3.Data$Analysis$PGStructs)]]
Proc.Trap_D3 <- PlotOnStages(Structure = "Circle",
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
                             SmoothPoints = 2,
                             MinCellPerNode = 2,
                             Title = 'Trapnell et al (D3)')


Staged.Trapnel.D3 <- StageWithPeaks(DataStruct = Trap_D3.Data,
                                    ProcStruct = Proc.Trap_D3,
                                    ComputeG0 = TRUE,
                                    FiltMax = 0,
                                    Thr = .8,
                                    QuantSel = .8,
                                    StageInfo = SelStageInfo.Human,
                                    MinNodes = 2,
                                    Mode = 1,
                                    G0Level = 1.1,
                                    Title = 'Trapnell et al D3')



















TB0 <- table(Staged.Trapnel.D0$CellStages_Ext)
TB1 <- table(Staged.Trapnel.D1$CellStages_Ext)
TB2 <- table(Staged.Trapnel.D2$CellStages_Ext)
TB3 <- table(Staged.Trapnel.D3$CellStages_Ext)

barplot(c(sum(TB0[grep('G0', names(TB0))])/sum(TB0),
          sum(TB1[grep('G0', names(TB1))])/sum(TB1),
          sum(TB2[grep('G0', names(TB2))])/sum(TB2),
          sum(TB3[grep('G0', names(TB3))])/sum(TB3)),
        names.arg = c("d0", "d1", "d2", "d3")
)


TB0 <- table(Staged.Trapnel.D0$CellStages)
TB1 <- table(Staged.Trapnel.D1$CellStages)
TB2 <- table(Staged.Trapnel.D2$CellStages)
TB3 <- table(Staged.Trapnel.D3$CellStages)

sum(TB0[grep('G0', names(TB0))])/sum(TB0)
sum(TB1[grep('G0', names(TB1))])/sum(TB1)
sum(TB2[grep('G0', names(TB2))])/sum(TB2)
sum(TB3[grep('G0', names(TB3))])/sum(TB3)


















D0Genes <- Output.Trap_D0$Genes[[length(Output.Trap_D0$Genes)]]
D1Genes <- Output.Trap_D1$Genes[[length(Output.Trap_D1$Genes)]]
D2Genes <- Output.Trap_D2$Genes[[length(Output.Trap_D2$Genes)]]
D3Genes <- Output.Trap_D3$Genes[[length(Output.Trap_D3$Genes)]]

AllGenes <- intersect(intersect(D0Genes, D1Genes), intersect(D2Genes, D3Genes))

length(D0Genes)/length(AllGenes)
length(D1Genes)/length(AllGenes)
length(D2Genes)/length(AllGenes)
length(D3Genes)/length(AllGenes)


PlotOnPseudotime(WorkStruct = InputList[[i]]$OrderedData,
                      Expression = InputList[[i]]$Expression,
                      Name = InputList[[i]]$Name,
                      gName = "Cdk4", SpanVal = SpanVect[i], CatOrder = NULL)



#
#
#
#
# TargetStruct <- Output.Trap_D0$PGStructs[[length(Output.Trap_D0$PGStructs)]]
# Proc.Exp.Trap <- PlotOnStages(Structure = "Circle",
#                               Categories = TargetStruct$Categories,
#                               nGenes = 2,
#                               TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
#                               PrinGraph = TargetStruct$PrinGraph,
#                               Net = TargetStruct$Net[[length(TargetStruct$Net)]],
#                               SelThr = .3,
#                               ComputeOverlaps = TRUE,
#                               ExpData = TargetStruct$FiltExp,
#                               RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
#                               PCACenter = TargetStruct$PCAData$center,
#                               PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
#                               OrderOnCat = FALSE,
#                               SmoothPoints = 2, MinCellPerNode = 2)
#
#
#
# pheatmap::pheatmap(Proc.Exp.Trap$NodesExp[rownames(Proc.Exp.Trap$NodesExp) %in% toupper(G1.Sel) |
#                                             rownames(Proc.Exp.Trap$NodesExp) %in% toupper(S.Sel) |
#                                             rownames(Proc.Exp.Trap$NodesExp) %in% toupper(G2.Sel)],
#                    cluster_cols = FALSE)
#
#
# plot(Proc.Exp.Trap$NodesExp["CDK4",])
#
#
#
#
# Output.Trap_D0$PGStructs[[5]]$PrinGraph
#
#
# ProcStruct <- Proc.Exp.Trap
# DataStruct <- list(Analysis = Output.Trap_D0,
#                    ExpMat = AllData)
#
#
#
# DataStruct = Data.Buet
# ProcStruct = Proc.Exp.Buet
# FiltMax = 0
# Thr = .7
# QuantSel = .5
# SinglePeack = TRUE
#
#
#
#
#   pheatmap::pheatmap(PlotMat[rownames(PlotMat) %in% G1.All,], show_colnames = FALSE,
#                      show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)
#
#   pheatmap::pheatmap(PlotMat[rownames(PlotMat) %in% S.All,], show_colnames = FALSE,
#                      show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)
#
#   pheatmap::pheatmap(PlotMat[rownames(PlotMat) %in% G2.All,], show_colnames = FALSE,
#                      show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)
#
#
#   # pheatmap::pheatmap(GeneExprMat.DF.Split.Mean.Bind[rownames(GeneExprMat.DF.Split.Mean.Bind) %in% G1.Sel,],
#   #                    cluster_rows = FALSE, cluster_cols = FALSE)
#   #
#   # pheatmap::pheatmap(GeneExprMat.DF.Split.Mean.Bind[rownames(GeneExprMat.DF.Split.Mean.Bind) %in% S.Sel,],
#   #                    cluster_rows = FALSE, cluster_cols = FALSE)
#   #
#   # pheatmap::pheatmap(GeneExprMat.DF.Split.Mean.Bind[rownames(GeneExprMat.DF.Split.Mean.Bind) %in% G2.Sel,],
#   #                    cluster_rows = FALSE, cluster_cols = FALSE)
#
#   InfoData <- rbind(
#     colSums(PlotMat[rownames(PlotMat) %in% G1.All,])/sum(rownames(PlotMat) %in% G1.All),
#     colSums(PlotMat[rownames(PlotMat) %in% S.All,])/sum(rownames(PlotMat) %in% S.All),
#     colSums(PlotMat[rownames(PlotMat) %in% G2.All,]/sum(rownames(PlotMat) %in% G2.All))
#   )
#
#   InfoData <- rbind(
#     colSums(PlotMat[rownames(PlotMat) %in% G1.Sel,]),
#     colSums(PlotMat[rownames(PlotMat) %in% S.Sel,]),
#     colSums(PlotMat[rownames(PlotMat) %in% G2.Sel,])
#   )
#
#   1*(InfoData > .04)
#
#   apply(InfoData, 2, which.max)
#
#   barplot(InfoData, beside = TRUE, las = 2)
#
#
#
#
#
#
#
#
#   StageGenes.Names <- lapply(StageGenes, function(x){
#     names(which(x))
#   })
#
#   barplot(unlist(lapply(StageGenes.Names, length)), las = 2, horiz = FALSE)
#
#   return(StageGenes.Names)
#
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# #
# # DataSet = AllData
# # StartSet = HumanGenes_GOCellCycle
# # Categories = factor(Categories)
# #
# # VarThr = .99
# # nNodes = 40
# # Log = TRUE
# # Filter = TRUE
# # OutThr = 5
# # PCAFilter = TRUE
# # OutThrPCA = 3
# # GraphType = "Circle"
# # PlanVarLimit = .85
# # PlanVarLimitIC = .9
# # MinBranDiff = 3
# # InitStructNodes = 20
# # ForceLasso = FALSE
# # EstProlif = "MeanPerc"
# # QuaThr = .5
# # NonG0Cell = NULL
# # DipPVThr = 1e-4
# # MinProlCells = 50
# # PCACenter = TRUE
# # PCAProjCenter = TRUE
# # PlotDebug = FALSE
# # PlotIntermediate = FALSE
# # AddGenePerc = 5
# # SelThr1 = .95
# # SelThr2 = .99
# # MadsThr =  1
#
#
#
#
