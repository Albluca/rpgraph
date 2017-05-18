# Load libraries ----------------------------------------------------------------

options(java.parameters = "-Xmx4g")

library(dplyr)
library(ggplot2)
library(rpgraph)
library(readr)
library(readxl)
library(GO.db)
library(biomaRt)


# Set geneset ----------------------------------------------------------------


FreemanData <- read_delim("~/Google Drive/Datasets/Gene List/Freeman.gmt", 
                          "\t", escape_double = FALSE, col_names = FALSE, 
                          trim_ws = TRUE)

Freeman_CC1 <- unlist(FreemanData[1,-c(1,2)], use.names = FALSE)
Freeman_CC2 <- unlist(FreemanData[2,-c(1,2)], use.names = FALSE)
Freeman_G1S_CC4 <- unlist(FreemanData[3,-c(1,2)], use.names = FALSE)
Freeman_G2M_CC6 <- unlist(FreemanData[4,-c(1,2)], use.names = FALSE)
Freeman_CC9 <- unlist(FreemanData[5,-c(1,2)], use.names = FALSE)
Freeman_G1S_CC4A <- unlist(FreemanData[6,-c(1,2)], use.names = FALSE)
Freeman_S_CC4B <- unlist(FreemanData[7,-c(1,2)], use.names = FALSE)
Freeman_G2_CC6A <- unlist(FreemanData[8,-c(1,2)], use.names = FALSE)
Freeman_M_CC6B <- unlist(FreemanData[9,-c(1,2)], use.names = FALSE)

Freeman_CC1_Mouse <- ConvertNames("human", "mouse", Freeman_CC1)
Freeman_CC2_Mouse <- ConvertNames("human", "mouse", Freeman_CC2)
Freeman_G1S_CC4_Mouse <- ConvertNames("human", "mouse", Freeman_G1S_CC4)
Freeman_G2M_CC6_Mouse <- ConvertNames("human", "mouse", Freeman_G2M_CC6)
Freeman_CC9_Mouse <- ConvertNames("human", "mouse", Freeman_CC9)
Freeman_G1S_CC4A_Mouse <- ConvertNames("human", "mouse", Freeman_G1S_CC4A)
Freeman_S_CC4B_Mouse <- ConvertNames("human", "mouse", Freeman_S_CC4B)
Freeman_G2_CC6A_Mouse <- ConvertNames("human", "mouse", Freeman_G2_CC6A)
Freeman_M_CC6B_Mouse <- ConvertNames("human", "mouse", Freeman_M_CC6B)

AllFreman_Mouse <- c(Freeman_CC1_Mouse, Freeman_CC2_Mouse, Freeman_CC9_Mouse,
                     Freeman_G1S_CC4_Mouse, Freeman_G1S_CC4A_Mouse, Freeman_G2_CC6A_Mouse,
                     Freeman_G2M_CC6_Mouse, Freeman_M_CC6B_Mouse, Freeman_S_CC4B_Mouse)

AllFreman_Mouse_CCP <- c(Freeman_G1S_CC4_Mouse, Freeman_G1S_CC4A_Mouse, Freeman_G2_CC6A_Mouse,
                         Freeman_G2M_CC6_Mouse, Freeman_M_CC6B_Mouse, Freeman_S_CC4B_Mouse)


AllTerms <- GOBPOFFSPRING[["GO:0007049"]]
MouseGenes_GOCellCycle <- getBM(attributes = c("external_gene_name"),
                                filters = "go",
                                values = AllTerms,
                                mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))

MouseGenes_GOCellCycle <- unlist(MouseGenes_GOCellCycle, use.names = FALSE)

AllWit_Mouse <- ConvertNames("human", "mouse", unlist(StageAssociation_Whit_Ext[-c(1,2)], use.names = FALSE))



# Define functions ----------------------------------------------------------------














# Buettner et al ----------------------------------------------------------------

# Load data

BaseDir <- "~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/"

Murine_ESC_G1_Data <- read_delim(paste(BaseDir, "G1_singlecells_counts.txt", sep= ''), "\t", escape_double = FALSE, trim_ws = TRUE)
Murine_ESC_S_Data <- read_delim(paste(BaseDir, "S_singlecells_counts.txt", sep= ''), "\t", escape_double = FALSE, trim_ws = TRUE)  
Murine_ESC_G2M_Data <- read_delim(paste(BaseDir, "G2M_singlecells_counts.txt", sep= ''), "\t", escape_double = FALSE, trim_ws = TRUE)

GeneToConsider <- which(!is.na(Murine_ESC_G1_Data$AssociatedGeneName))

Murine_ESC_All <- cbind(Murine_ESC_G1_Data[GeneToConsider, ],
                        Murine_ESC_S_Data[GeneToConsider, -c(1:4)],
                        Murine_ESC_G2M_Data[GeneToConsider, -c(1:4)])

GNames <- Murine_ESC_All$AssociatedGeneName

GNames[duplicated(GNames)] <- paste(GNames[duplicated(GNames)], 1:sum(duplicated(GNames)), sep ="_R")

Murine_ESC_All <- Murine_ESC_All[,-c(1:4)]
rownames(Murine_ESC_All) <- GNames

Murine_ESC_Stages <- c(rep(x = "G1", ncol(Murine_ESC_G1_Data)-4), rep(x = "S", ncol(Murine_ESC_S_Data)-4), rep(x = "G2M", ncol(Murine_ESC_G2M_Data)-4))
names(Murine_ESC_Stages) <- colnames(Murine_ESC_All)


Murine_ESC_Stages <- factor(Murine_ESC_Stages, levels = c("G1", "S", "G2M"))

GeneExprMat <- Murine_ESC_All
StartSet <- MouseGenes_GOCellCycle
Categories <- Murine_ESC_Stages


# BuettInfo <- PlotsAndDistill(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .5,
#                                    IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9,
#                                    PCACenter = FALSE, PlotDebug = TRUE, CompareSet = list("Whit" = AllWit_Mouse, "Free" = AllFreman_Mouse))
# 
# BuettProc <- PlotOnStages(Structure = "Circle", TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]],
#                                Categories = BuettInfo$FinalStruct$Categories, nGenes = 10,
#                                PrinGraph = BuettInfo$FinalStruct$PrinGraph,
#                                Net = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
#                                SelThr = .35, ComputeOverlaps = TRUE, ExpData = BuettInfo$FinalStruct$FiltExp,
#                                RotatioMatrix = BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims],
#                                PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]])
# 
# 
# 
# ColCat <- c("red", "blue", "green")
# names(ColCat) <- c("G1", "S", "G2M")
# 
# 
# plotPieNet(Results =  BuettInfo$FinalStruct$IntGrahs[[length(BuettInfo$FinalStruct$IntGrahs)]],
#            Data = BuettInfo$FinalStruct$Data, NodeSizeMult = 4,
#            Categories = BuettInfo$FinalStruct$Categories, PlotNet = TRUE,
#            Graph = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
#            TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]],
#            LayOut = 'circle', Main = "Pincipal graph", ColCat = ColCat)
# 
# legend(x = "center", legend = names(ColCat), fill = ColCat)

# BuettProc <- list()
# BuettInfoList <- list()
# 
# 
# for(i in 1:10){
#   
#   BuettInfo <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .7,
#                                      IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9)
#   
#   BuettInfoList[[i]] <- BuettInfo
#   
#   
#   BuettProc[[i]] <- PlotOnStages(Structure = "Circle",
#                                  TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]],
#                                  Categories = BuettInfo$FinalStruct$Categories, nGenes = 2,
#                                  PrinGraph = BuettInfo$FinalStruct$PrinGraph,
#                                  Net = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
#                                  SelThr = .35, ComputeOverlaps = TRUE, 
#                                  ExpData = BuettInfo$FinalStruct$FiltExp,
#                                  RotatioMatrix = BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims],
#                                  PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]])
#   
# }
# 
# 
# AllCellsPT <- lapply(BuettProc, "[[", "CellsPT")
# 
# ReordCells <- sapply(AllCellsPT, function(x) {
#   tData <- rep(NA, ncol(GeneExprMat))
#   names(tData) <- colnames(GeneExprMat)
#   tData[names(x)] <- x
#   tData
# })
# 
# boxplot(t(apply(ReordCells, 2, rank, ties.method = "min"))[,order(apply(ReordCells, 1, median))],
#         xlab = "Cells", ylab = "Position", xaxt='n')
# 
# 
# plot(apply(apply(ReordCells, 2, rank, ties.method = "min"), 1, median)[order(apply(ReordCells, 1, median))],
#      ylab = "Median order", xlab = "Cells")
# 
# 
# hist(apply(apply(ReordCells, 2, rank, ties.method = "min"), 1, var))
# 
# 
# 
# 
# 
# 
# AllGenes <- lapply(BuettProc, function(x){
#   rownames(x$NodesExp)
# })
# 
# AllGenesNales <- unique(unlist(AllGenes, use.names = FALSE))
# 
# GeneOrder <- sapply(AllGenes, function(x) {
#   match(AllGenesNales, x)
# })
# 
# boxplot(t(GeneOrder[order(apply(GeneOrder, 1, median, na.rm=TRUE)),]))



























# for(j in 1:10){
#   
#   Gname <- sample(rownames(BuettProc[[1]]$NodesExp), 1)
#   # Gname <- rownames(BuettProc[[1]]$NodesExp)[1]
#   
#   tData <- NULL
#   
#   for(i in 1:length(BuettProc)){
#     
#     if(!any(rownames(BuettProc[[i]]$NodesExp) == Gname))
#       next()
#     
#     tData <- rbind(tData, cbind(BuettProc[[i]]$NodesPT/sum(BuettInfoList[[i]]$FinalStruct$ProjPoints[[length(BuettInfoList[[i]]$FinalStruct$ProjPoints)]]$EdgeLength),
#                                 BuettProc[[i]]$NodesExp[Gname, ],
#                                 rep(i, length(BuettProc[[i]]$NodesPT))))
#   }
#   
#   colnames(tData) <- c("PT", "Exp", "Rep")
#   tDF <- data.frame(tData)
#   tDF$Rep <- as.character(tDF$Rep)
#   
#   p <- ggplot(data = tDF, mapping = aes(x = PT, y = Exp)) + geom_smooth(span=.25) +
#     geom_point(mapping = aes(color = Rep)) + labs(title = Gname)
#   
#   print(p)
#   
#   
# }
# 
# 
# RectData <- lapply(as.list(1:length(BuettProc)), function(i){
#   tLen <- sum(BuettInfoList[[i]]$FinalStruct$ProjPoints[[length(BuettInfoList[[i]]$FinalStruct$ProjPoints)]]$EdgeLength)
#   cbind(BuettProc[[i]]$RecCoord$Min/tLen,
#         BuettProc[[i]]$RecCoord$Med/tLen,
#         BuettProc[[i]]$RecCoord$Max/tLen,
#         as.character(BuettProc[[i]]$RecCoord$Stage),
#         rep(i, length(BuettProc[[i]]$RecCoord$Stage)))}
# )
#   
#   
# RectDF <- NULL
# 
# for(i in 1:length(RectData)){
#   RectDF <- rbind(RectDF, RectData[[i]])
# }
# 
# colnames(RectDF) <- c("Min", "Med", "Max", "Stage", "Low")
# 
# RectDF.Final <- data.frame(Min=as.numeric(RectDF[,1]),
#                            Med=as.numeric(RectDF[,2]),
#                            Max=as.numeric(RectDF[,3]),
#                            Low = as.numeric(RectDF[,5]),
#                            High = as.numeric(RectDF[,5]) + 1,
#                            Stage = RectDF[,4])
# 
# 
# 
# 
# 
# 
# # Gname <- sample(rownames(BuettProc[[1]]$NodesExp), 1)
# Gname <- rownames(BuettProc[[1]]$NodesExp)[1]
# 
# tData <- NULL
# 
# for(i in 1:length(BuettProc)){
#   
#   if(!any(rownames(BuettProc[[i]]$NodesExp) == Gname))
#     next()
#   
#   tData <- rbind(tData, cbind(BuettProc[[i]]$NodesPT/sum(BuettInfoList[[i]]$FinalStruct$ProjPoints[[length(BuettInfoList[[i]]$FinalStruct$ProjPoints)]]$EdgeLength),
#                               BuettProc[[i]]$NodesExp[Gname, ],
#                               rep(i, length(BuettProc[[i]]$NodesPT))))
# }
# 
# colnames(tData) <- c("PT", "Exp", "Rep")
# tDF <- data.frame(tData)
# tDF$Rep <- as.character(tDF$Rep)
# 
# 
# Loc.RectDF.Final <- RectDF.Final
# Loc.RectDF.Final$Low <- Loc.RectDF.Final$Low - 1
# Loc.RectDF.Final$High <- Loc.RectDF.Final$High - 1
# 
# MinExp <- min(tDF$Exp)
# MaxExp <- max(tDF$Exp)
# RangeExp <- MaxExp - MinExp
# StepCat <- RangeExp/max(Loc.RectDF.Final$High)
# 
# Loc.RectDF.Final$Low <- Loc.RectDF.Final$Low*StepCat + MinExp
# Loc.RectDF.Final$High <- Loc.RectDF.Final$High*StepCat + MinExp
# 
# p <- ggplot(data = tDF, mapping = aes(x = PT, y = Exp)) +
#   geom_rect(data = Loc.RectDF.Final, mapping = aes(xmin = Min, xmax = Max, ymin = Low, ymax = High, fill = Stage),
#             inherit.aes = FALSE, alpha = .3) + geom_smooth(span=.25) +
#   geom_point(mapping = aes(color = Rep)) + labs(title = Gname)
# 
# print(p)



















# RankVar <- apply(apply(ReordCells, 2, rank, ties.method = "min"), 1, var)
# RankVar[which(RankVar > mean(RankVar))]
# 
# KM <- kmeans(RankVar, centers = range(RankVar))
# if(KM$centers[1] < KM$centers[2]){
#   ToUse <- which(KM$cluster == 1)
# } else {
#   ToUse <- which(KM$cluster == 2)
# }













# AllGenes <- lapply(BuettInfoList, function(x){
#   colnames(x$FinalStruct$FiltExp)
# })
# 
# AllGenesNales <- unique(unlist(AllGenes, use.names = FALSE))
# 
# GeneOrder <- sapply(AllGenes, function(x) {
#   match(AllGenesNales, x)
# })
# 
# 
# AllGenesNales[which(rowSums(is.na(GeneOrder)) == 0)]
# 

# AllFiltGenes <- rownames(GeneExprMat)
# 
# lapply(BuettInfoList, function(x){
#   AllFiltGenes <<- intersect(AllFiltGenes, x$Genes)
# })



# FilteredFinalList <- list()
# BuettProcList <- list()
# 
# for(i in 1:10){
#   FilteredFinal <- ProjectAndCompute(DataSet = GeneExprMat,
#                                      GeneSet = AllGenesNales[which(rowSums(is.na(GeneOrder)) == 0)],
#                                      Categories = Categories, nNodes = 40,
#                                      VarThr = .99, GraphType = 'Circle', PlanVarLimit = .85,
#                                      PlanVarLimitIC = .9, ForceLasso = FALSE, InitStructNodes = 20,
#                                      MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 20, DipPVThr = 1e-4,
#                                      PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE, PCAFilter = TRUE,
#                                      PlotDebug = FALSE)
#   
#   FilteredFinalList[[i]] <- FilteredFinal
#   
#   ProjectOnPrincipalGraph(Nodes = FilteredFinal$PrinGraph$Nodes, Edges = FilteredFinal$PrinGraph$Edges, Points = FilteredFinal$Data,
#                           UsedPoints = NULL, Categories = FilteredFinal$Categories, Title=paste("(Filtered)"),
#                           PCACenter = TRUE, ShowFitted = FALSE)
#   
#   BuettProc <- PlotOnStages(Structure = "Circle", TaxonList = FilteredFinal$TaxonList[[length(FilteredFinal$TaxonList)]],
#                             Categories = FilteredFinal$Categories, nGenes = 2,
#                             PrinGraph = FilteredFinal$PrinGraph,
#                             Net = FilteredFinal$Net[[length(FilteredFinal$Net)]],
#                             SelThr = .35, ComputeOverlaps = TRUE, ExpData = FilteredFinal$FiltExp,
#                             RotatioMatrix = FilteredFinal$PCAData$rotation[,1:FilteredFinal$nDims],
#                             PointProjections = FilteredFinal$ProjPoints[[length(FilteredFinal$ProjPoints)]])
#   
#   BuettProcList[[i]] <- BuettProc
# }


















# FilteredFinal <- ProjectAndCompute(DataSet = GeneExprMat,
#                                    GeneSet = AllGenesNales[which(rowSums(is.na(GeneOrder)) == 0)],
#                                    Categories = Categories, nNodes = 40,
#                                    VarThr = .99, GraphType = 'Circle', PlanVarLimit = .85,
#                                    PlanVarLimitIC = .9, ForceLasso = FALSE, InitStructNodes = 20,
#                                    MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 20, DipPVThr = 1e-4,
#                                    PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE, PCAFilter = TRUE,
#                                    PlotDebug = FALSE)
# 
# ProjectOnPrincipalGraph(Nodes = FilteredFinal$PrinGraph$Nodes, Edges = FilteredFinal$PrinGraph$Edges, Points = FilteredFinal$Data,
#                         UsedPoints = NULL, Categories = FilteredFinal$Categories, Title=paste("(Filtered)"),
#                         PCACenter = TRUE, ShowFitted = FALSE)
# 
# BuettProc <- PlotOnStages(Structure = "Circle", TaxonList = FilteredFinal$TaxonList[[length(FilteredFinal$TaxonList)]],
#                           Categories = FilteredFinal$Categories, nGenes = 2,
#                           PrinGraph = FilteredFinal$PrinGraph,
#                           Net = FilteredFinal$Net[[length(FilteredFinal$Net)]],
#                           SelThr = .35, ComputeOverlaps = TRUE, ExpData = FilteredFinal$FiltExp,
#                           RotatioMatrix = FilteredFinal$PCAData$rotation[,1:FilteredFinal$nDims],
#                           PointProjections = FilteredFinal$ProjPoints[[length(FilteredFinal$ProjPoints)]])















# 
# # Exapand genes
# 
# # WORK HERE!!!!
# 
# 
# 
# Topology = 'Circle'
# DistillThr = .6
# IgnoreTail = TRUE
# Log = TRUE
# StartQuant = .25
# Title = "Buettner et al"
# PlanVarLimit = .85
# PlanVarLimitIC = .9
# KeepOriginal = TRUE
# PCACenter = FALSE
# PlotDebug = FALSE
# Mode = "VarPC"
# ExtMode = 3
# MaxRounds = 15
# 








GeneExprMat = Murine_ESC_All
StartSet = MouseGenes_GOCellCycle
Categories = Murine_ESC_Stages
Topology = 'Circle'
DistillThr = .5
IgnoreTail = FALSE
Log = TRUE
StartQuant = .5
Title = "Buettner et al. (2015) "
PlanVarLimit = .85
PlanVarLimitIC = .9
KeepOriginal = TRUE
PCACenter = FALSE
PlotDebug = FALSE
QuaThr = .5
Mode = "VarPC"
ExtMode = 2
StopCrit = .95
StopMode = c(1, 2)
ExpQuant = .01
Parallel = TRUE
nCores = 3
ClusType = "FORK"
MinProlCells = 25
SelThr = .3
KeepOrgProlifiration = TRUE
Filter = TRUE

OutThrPCA = 3
EstProlif = "MeanPerc"
PCAFilter = TRUE
OutThr = 3
Filter = TRUE
InitStructNodes = 20
nNodes = 40
PCAProjCenter = TRUE


ExpData = GeneExprMat
StartGeneSet = StartSet
VarThr = .99
GraphType = 'Circle'
PlotIntermediate = FALSE


# BaseAnalysis = ConsideredStruct
# FullExpression = ExpData
# Topo = 'Circle'



ExpDataset.Buett <- ProjectAndExpand(GeneExprMat = Murine_ESC_All, StartSet = MouseGenes_GOCellCycle, Categories = Murine_ESC_Stages,
                               Topology = 'Circle', DistillThr = .5, IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Buettner et al. (2015) ",
                               PlanVarLimit = .85, PlanVarLimitIC = .9, KeepOriginal = TRUE, PCACenter = FALSE, PlotDebug = FALSE, QuaThr = .5,
                               Mode = "VarPC", ExtMode = 2, StopCrit = .95, StopMode = c(1, 2), ExpQuant = .05, 
                               Parallel = TRUE, nCores = 3, ClusType = "FORK", MinProlCells = 25, SelThr = .35,
                               OrderOnCat = TRUE, KeepOrgProlifiration = TRUE)

# 
# ExpDataset.Kowa <- ProjectAndExpand(GeneExprMat = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
#                                     Topology = 'Circle', DistillThr = .5, IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al. (2015)",
#                                     PlanVarLimit = .85, PlanVarLimitIC = .9, KeepOriginal = TRUE, PCACenter = FALSE, PlotDebug = FALSE, QuaThr = .5,
#                                     Mode = "VarPC", ExtMode = 2, StopCrit = .98, StopMode = c(1, 2), ExpQuant = .05, EstProlif = "MeanPerc",
#                                     Parallel = TRUE, nCores = 3, ClusType = "FORK", MinProlCells = 50, SelThr = .35,
#                                     OrderOnCat = TRUE, KeepOrgProlifiration = TRUE)




FullExpData.Buett <- log10(Murine_ESC_All + 1)
FullCat.Buett <- Murine_ESC_Stages

save(ExpDataset.Buett, FullExpData.Buett, FullCat.Buett,
     ExpDataset.Kowa, FullExpData.Kowa, FullCat.Kowa,
     ExpDataset.Sasa, FullExpData.Sasa, FullCat.Sasa,
     file = "/Users/newmac-luca/Desktop/RecentData_U.RData")




BuettGenes <- ExpDataset.Buett$GeneList[[length(ExpDataset.Buett$GeneList)]]

FinalStructure <- ExpDataset.Buett$ListProc[[length(ExpDataset.Buett$ListProc)]]

ReOrd <- match(colnames(FinalStructure$CellExp), names(FinalStructure$CellsPT))

# plot(FinalStructure$CellsPT[ReOrd], FinalStructure$CellExp["Mcm5",])
# points(FinalStructure$NodesPT, FinalStructure$NodesExp["Mcm5",], col='red', type = "b", lwd=3)



Genes.Buett <- sort(unique(c(rownames(FinalStructure$NodesExp), rownames(FinalStructure$NodesExp))))

# gName <- "Cdk4"


p <- ggplot2::ggplot(data = data.frame(x=FinalStructure$NodesPT,
                                       y=FinalStructure$NodesExp[gName,]),
                     mapping = ggplot2::aes(x = x, y = y, color="PC")) +
  ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(gName, " - Buettner et al.", sep ="")) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))

p <- p + ggplot2::geom_rect(data = FinalStructure$RecCoord,
                            mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                            ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
  ggplot2::geom_point(data = data.frame(x=FinalStructure$CellsPT[ReOrd],
                                        y=FinalStructure$CellExp[gName,]),
                      mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
  ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))

print(p)


















GenesOnNodes <- ExpDataset$ListProc[[length(ExpDataset$ListProc)]]$NodesExp
Q95 <- apply(GenesOnNodes, 1, quantile, .95)
# View((GenesOnNodes > Q95)*1)


StagesList <- ExpDataset$ListProc[[length(ExpDataset$ListProc)]]$StageOnNodes

for(i in 1:length(StagesList)){
  if(length(StagesList[[i]]) > 0){
    StagesList[[i]] <- sort(unique(c(StagesList[[i]],
                                     min(StagesList[[i]]):max(StagesList[[i]]))))
  }
}

MaxNode <- apply(GenesOnNodes, 1, which.max)
barplot(table(MaxNode), log='y')


StagedGenes <- lapply(StagesList, function(x){
  MaxNode[MaxNode %in% x]
})

barplot(unlist(lapply(StagedGenes, length)))

StageAssociation_Whit$Stages
StageAssociation_Whit.Mouse.S1 <- ConvertNames(SourceOrganism = "human", TargetOrganism = "mouse", StageAssociation_Whit$S1_U)
StageAssociation_Whit.Mouse.S2 <- ConvertNames(SourceOrganism = "human", TargetOrganism = "mouse", StageAssociation_Whit$S2_U)
StageAssociation_Whit.Mouse.S3 <- ConvertNames(SourceOrganism = "human", TargetOrganism = "mouse", StageAssociation_Whit$S3_U)
StageAssociation_Whit.Mouse.S4 <- ConvertNames(SourceOrganism = "human", TargetOrganism = "mouse", StageAssociation_Whit$S4_U)


StageAssociation




intersect(unlist(lapply(StagedGenes, names)),
          c(StageAssociation_Whit.Mouse.S1, StageAssociation_Whit.Mouse.S2,
            StageAssociation_Whit.Mouse.S3, StageAssociation_Whit.Mouse.S4))



intersect(c(names(StagedGenes$G1), names(StagedGenes$`G1+S`)),
          StageAssociation_Whit.Mouse.S1)
intersect(c(names(StagedGenes$G1), names(StagedGenes$`G1+S`)),
          StageAssociation_Whit.Mouse.S2)
intersect(c(names(StagedGenes$G1), names(StagedGenes$`G1+S`)),
          StageAssociation_Whit.Mouse.S3)
intersect(c(names(StagedGenes$G1), names(StagedGenes$`G1+S`)),
          StageAssociation_Whit.Mouse.S4)


intersect(names(StagedGenes$S), StageAssociation_Whit.Mouse.S1)
intersect(names(StagedGenes$S), StageAssociation_Whit.Mouse.S2)
intersect(names(StagedGenes$S), StageAssociation_Whit.Mouse.S3)
intersect(names(StagedGenes$S), StageAssociation_Whit.Mouse.S4)




intersect(names(StagedGenes$G2M), StageAssociation_Whit.Mouse.S1)
intersect(names(StagedGenes$G2M), StageAssociation_Whit.Mouse.S2)
intersect(names(StagedGenes$G2M), StageAssociation_Whit.Mouse.S3)
intersect(names(StagedGenes$G2M), StageAssociation_Whit.Mouse.S4)





# TO DO: Look at peacks










# BuettInfo.Exp.75 <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .75))), MouseGenes_GOCellCycle)),
#                                           Categories = Categories, Topology = 'Circle', DistillThr = .5,
#                                           IgnoreTail = TRUE, Log = TRUE, StartQuant = .75, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9,
#                                           PCACenter = FALSE, PlotDebug = FALSE)
# 
# BuettInfo.Exp.25 <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .25))), MouseGenes_GOCellCycle)),
#                                           Categories = Categories, Topology = 'Circle', DistillThr = .5,
#                                           IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9,
#                                           PCACenter = FALSE, PlotDebug = FALSE)
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
# BuettProc.Exp.25 <- PlotOnStages(Structure = "Circle", TaxonList = BuettInfo.Exp.25$FinalStruct$TaxonList[[length(BuettInfo.Exp.25$FinalStruct$TaxonList)]],
#                                 Categories = BuettInfo.Exp.25$FinalStruct$Categories, nGenes = 2,
#                                 PrinGraph = BuettInfo.Exp.25$FinalStruct$PrinGraph,
#                                 Net = BuettInfo.Exp.25$FinalStruct$Net[[length(BuettInfo.Exp.25$FinalStruct$Net)]],
#                                 SelThr = .35, ComputeOverlaps = TRUE, ExpData = BuettInfo.Exp.25$FinalStruct$FiltExp,
#                                 RotatioMatrix = BuettInfo.Exp.25$FinalStruct$PCAData$rotation[,1:BuettInfo.Exp.25$FinalStruct$nDims],
#                                 PointProjections = BuettInfo.Exp.25$FinalStruct$ProjPoints[[length(BuettInfo.Exp.25$FinalStruct$ProjPoints)]])
# 
# 
# BuettProc.Exp.75 <- PlotOnStages(Structure = "Circle", TaxonList = BuettInfo.Exp.75$FinalStruct$TaxonList[[length(BuettInfo.Exp.75$FinalStruct$TaxonList)]],
#                                  Categories = BuettInfo.Exp.75$FinalStruct$Categories, nGenes = 2,
#                                  PrinGraph = BuettInfo.Exp.75$FinalStruct$PrinGraph,
#                                  Net = BuettInfo.Exp.75$FinalStruct$Net[[length(BuettInfo.Exp.75$FinalStruct$Net)]],
#                                  SelThr = .35, ComputeOverlaps = TRUE, ExpData = BuettInfo.Exp.75$FinalStruct$FiltExp,
#                                  RotatioMatrix = BuettInfo.Exp.75$FinalStruct$PCAData$rotation[,1:BuettInfo.Exp.75$FinalStruct$nDims],
#                                  PointProjections = BuettInfo.Exp.75$FinalStruct$ProjPoints[[length(BuettInfo.Exp.75$FinalStruct$ProjPoints)]])
# 
# 
# 
# 
# 
# 
# 
# barplot(rbind(c(length(BuettInfo.Exp$Genes), length(BuettInfo.Exp.25$Genes), length(BuettInfo.Exp.5$Genes), length(BuettInfo.Exp.75$Genes)),
#               c(length(MouseGenes_GOCellCycle), length(unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .25))), MouseGenes_GOCellCycle))),
#                 length(unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .5))), MouseGenes_GOCellCycle))),
#                 length(unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .75))), MouseGenes_GOCellCycle)))))
# , beside = TRUE)
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
# sum(BuettInfo$Genes %in% BuettInfo.Exp$Genes)/length(BuettInfo$Genes)
# length(BuettInfo.Exp$Genes)/length(BuettInfo$Genes)
# 
# 
# 
# 
# CellToPlot <- intersect(names(BuettProc.Exp$CellsPT), names(BuettProc$CellsPT))
# 
# 
# plot(rank(BuettProc.Exp$CellsPT[CellToPlot], ties.method = "min"), rank(BuettProc$CellsPT[CellToPlot], ties.method = "min"))
# 
# 
# cor(BuettProc.Exp$CellsPT[CellToPlot], BuettProc$CellsPT[CellToPlot])
# 
# 
# cor(rank(BuettProc.Exp$CellsPT[CellToPlot], ties.method = "min"), rank(BuettProc$CellsPT[CellToPlot], ties.method = "min"))
# 
# 
# GenesByStage <- lapply(BuettProc$StageOnNodes, function(x){
#   BuettProc$NodesExp[,x]
# })
# 
# StageOnNodes <- rep(NA, ncol(BuettProc$NodesExp))
# lapply(as.list(1:length(BuettProc$StageOnNodes)),function(i){
#   StageOnNodes[BuettProc$StageOnNodes[[i]]] <<- names(BuettProc$StageOnNodes)[i]
# })
# StageOnNodes <- factor(StageOnNodes, levels = names(BuettProc$StageOnNodes))
# 
# # boxplot(BuettProc$NodesExp[1,] ~ StageOnNodes)
# 
# AllPV <- apply(BuettProc$NodesExp, 1, function(x){
#   AOV <- aov(x ~ StageOnNodes)
#   summary(AOV)[[1]][1,5]
# })
# 
# 
# TLGenes <- apply(BuettProc$NodesExp[which(AllPV < 1e-3),], 1, function(x){
#   AGG <- aggregate(x, by = list(StageOnNodes), median)
#   TS <- AGG[which.max(AGG[,2]),1]
#   LS <- AGG[which.min(AGG[,2]),1]
# 
#   c(as.character(TS), wilcox.test(x ~ StageOnNodes == TS, alternative = "less")$p.value,
#     as.character(LS), wilcox.test(x ~ StageOnNodes == LS, alternative = "greater")$p.value)
# 
# })
# 
# TLGenes <- t(TLGenes)
# 
# table(TLGenes[TLGenes[,4] < 1e-3, 3])
# table(TLGenes[TLGenes[,2] < 1e-3, 1])
# 
# # 
# # 
# # 
# # AllPV <- apply(BuettProc$NodesExp, 1, function(x){
# #   AOV <- aov(x ~ StageOnNodes)
# #   PV <- summary(AOV)[[1]][1,5]
# #   THT <- TukeyHSD(AOV)
# #   SigDiff <- which(THT$StageOnNodes[,"p adj"] < 1e-3)
# #   
# #   RestTHT <- THT$StageOnNodes[SigDiff,]
# #   SplGr <- strsplit(rownames(RestTHT), "-")
# #   Gr1 <- unlist(lapply(SplGr, "[[", 1))
# #   Gr2 <- unlist(lapply(SplGr, "[[", 2))
# #   DifSgn <- sign(RestTHT[, "diff"])
# # })
# # 
# # SelGenes <- which(AllPV < 1e-3)
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
# # 
