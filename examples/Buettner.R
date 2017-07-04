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
                                filters = "go_id",
                                values = AllTerms,
                                mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))

MouseGenes_GOCellCycle <- unlist(MouseGenes_GOCellCycle, use.names = FALSE)

AllWit_Mouse <- ConvertNames("human", "mouse", unlist(StageAssociation_Whit_Ext[-c(1,2)], use.names = FALSE))



# Define functions ----------------------------------------------------------------





PlotsAndDistill <- function(GeneExprMat, StartSet, Categories, Topology = "Circle", IgnoreTail = FALSE, PlanVarLimit = .9, PlanVarLimitIC = 92,
                                  DistillThr = .7, Log = TRUE, StartQuant = .5, OutThr = 3, OutThrPCA = 3, Title = '', MinProlCells = 20, PCACenter = FALSE,
                                  PlotDebug = FALSE, Mode = "VarPC", ExtMode = 2, nNodes = 40, InitStructNodes = 20, DipPVThr = 1e-4, PCAFilter = TRUE,
                                  PCAProjCenter = TRUE, StopCrit = .95, Filter = TRUE, CompareSet = list()) {
  
  
  # Perform the initial analysis
  
  PrGraph.Initial <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = StartSet, OutThr = OutThr, nNodes = nNodes, 
                                       VarThr = .99, Categories = Categories, GraphType = Topology, PlanVarLimit = PlanVarLimit,
                                       PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, InitStructNodes = InitStructNodes,
                                       MinBranDiff = 2, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                       PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = FALSE,
                                       PCAFilter = PCAFilter, OutThrPCA = OutThrPCA, PlotDebug = PlotDebug)
  
  # Distill genes
  
  DistilledGenes <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = Mode, DistillThr = DistillThr, FastReduce = FALSE,
                                    ExtMode = ExtMode, StopCrit = StopCrit, OutThr = OutThr, nNodes = nNodes, VarThr = .99, Categories = Categories,
                                    GraphType = Topology, PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE,
                                    InitStructNodes = InitStructNodes, MinBranDiff = 2, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                    PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = FALSE, PlotDebug = PlotDebug,
                                    IgnoreTail = IgnoreTail, StartQuant = StartQuant, PCAFilter = PCAFilter, OutThrPCA = OutThrPCA)
  
  # Perform the analysis on the final structure
  
  PrGraph.Final <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = DistilledGenes, OutThr = OutThr, nNodes = nNodes, 
                                     VarThr = .99, Categories = Categories, GraphType = Topology, PlanVarLimit = PlanVarLimit,
                                     PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, InitStructNodes = InitStructNodes,
                                     MinBranDiff = 2, Log = Log, Filter = Filter, MinProlCells = MinProlCells, DipPVThr = DipPVThr,
                                     PCACenter = PCACenter, PCAProjCenter = PCAProjCenter, PlotIntermediate = FALSE, PCAFilter = PCAFilter,
                                     OutThrPCA = OutThrPCA, PlotDebug = PlotDebug)
  
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
  
  if(length(CIListEnd) > 1){
    
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


BuettInfo <- PlotsAndDistill(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .5,
                                   IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9,
                                   PCACenter = FALSE, PlotDebug = TRUE, CompareSet = list("Whit" = AllWit_Mouse, "Free" = AllFreman_Mouse))

BuettProc <- PlotOnStages(Structure = "Circle", TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]],
                               Categories = BuettInfo$FinalStruct$Categories, nGenes = 10,
                               PrinGraph = BuettInfo$FinalStruct$PrinGraph,
                               Net = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
                               SelThr = .35, ComputeOverlaps = TRUE, ExpData = BuettInfo$FinalStruct$FiltExp,
                               RotatioMatrix = BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims],
                               PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]])



ColCat <- c("red", "blue", "green")
names(ColCat) <- c("G1", "S", "G2M")


plotPieNet(Results =  BuettInfo$FinalStruct$IntGrahs[[length(BuettInfo$FinalStruct$IntGrahs)]],
           Data = BuettInfo$FinalStruct$Data, NodeSizeMult = 4,
           Categories = BuettInfo$FinalStruct$Categories, PlotNet = TRUE,
           Graph = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
           TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]],
           LayOut = 'circle', Main = "Pincipal graph", ColCat = ColCat)

legend(x = "center", legend = names(ColCat), fill = ColCat)


# Structure = "Circle"
# TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]]
# Categories = BuettInfo$FinalStruct$Categories
# nGenes = 10
# PrinGraph = BuettInfo$FinalStruct$PrinGraph
# Net = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]]
# SelThr = .35
# ComputeOverlaps = TRUE
# ExpData = BuettInfo$FinalStruct$FiltExp
# RotatioMatrix = BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims]
# PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]]












BuettProc <- list()
BuettInfoList <- list()


for(i in 1:10){
  
  BuettInfo <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .7,
                                     IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9)
  
  BuettInfoList[[i]] <- BuettInfo
  
  
  BuettProc[[i]] <- PlotOnStages(Structure = "Circle",
                                 TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]],
                                 Categories = BuettInfo$FinalStruct$Categories, nGenes = 2,
                                 PrinGraph = BuettInfo$FinalStruct$PrinGraph,
                                 Net = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
                                 SelThr = .35, ComputeOverlaps = TRUE, 
                                 ExpData = BuettInfo$FinalStruct$FiltExp,
                                 RotatioMatrix = BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims],
                                 PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]])
  
}


AllCellsPT <- lapply(BuettProc, "[[", "CellsPT")

ReordCells <- sapply(AllCellsPT, function(x) {
  tData <- rep(NA, ncol(GeneExprMat))
  names(tData) <- colnames(GeneExprMat)
  tData[names(x)] <- x
  tData
})

boxplot(t(apply(ReordCells, 2, rank, ties.method = "min"))[,order(apply(ReordCells, 1, median))],
        xlab = "Cells", ylab = "Position", xaxt='n')


plot(apply(apply(ReordCells, 2, rank, ties.method = "min"), 1, median)[order(apply(ReordCells, 1, median))],
     ylab = "Median order", xlab = "Cells")


hist(apply(apply(ReordCells, 2, rank, ties.method = "min"), 1, var))






AllGenes <- lapply(BuettProc, function(x){
  rownames(x$NodesExp)
})

AllGenesNales <- unique(unlist(AllGenes, use.names = FALSE))

GeneOrder <- sapply(AllGenes, function(x) {
  match(AllGenesNales, x)
})

boxplot(t(GeneOrder[order(apply(GeneOrder, 1, median, na.rm=TRUE)),]))



























for(j in 1:10){
  
  Gname <- sample(rownames(BuettProc[[1]]$NodesExp), 1)
  # Gname <- rownames(BuettProc[[1]]$NodesExp)[1]
  
  tData <- NULL
  
  for(i in 1:length(BuettProc)){
    
    if(!any(rownames(BuettProc[[i]]$NodesExp) == Gname))
      next()
    
    tData <- rbind(tData, cbind(BuettProc[[i]]$NodesPT/sum(BuettInfoList[[i]]$FinalStruct$ProjPoints[[length(BuettInfoList[[i]]$FinalStruct$ProjPoints)]]$EdgeLength),
                                BuettProc[[i]]$NodesExp[Gname, ],
                                rep(i, length(BuettProc[[i]]$NodesPT))))
  }
  
  colnames(tData) <- c("PT", "Exp", "Rep")
  tDF <- data.frame(tData)
  tDF$Rep <- as.character(tDF$Rep)
  
  p <- ggplot(data = tDF, mapping = aes(x = PT, y = Exp)) + geom_smooth(span=.25) +
    geom_point(mapping = aes(color = Rep)) + labs(title = Gname)
  
  print(p)
  
  
}


RectData <- lapply(as.list(1:length(BuettProc)), function(i){
  tLen <- sum(BuettInfoList[[i]]$FinalStruct$ProjPoints[[length(BuettInfoList[[i]]$FinalStruct$ProjPoints)]]$EdgeLength)
  cbind(BuettProc[[i]]$RecCoord$Min/tLen,
        BuettProc[[i]]$RecCoord$Med/tLen,
        BuettProc[[i]]$RecCoord$Max/tLen,
        as.character(BuettProc[[i]]$RecCoord$Stage),
        rep(i, length(BuettProc[[i]]$RecCoord$Stage)))}
)
  
  
RectDF <- NULL

for(i in 1:length(RectData)){
  RectDF <- rbind(RectDF, RectData[[i]])
}

colnames(RectDF) <- c("Min", "Med", "Max", "Stage", "Low")

RectDF.Final <- data.frame(Min=as.numeric(RectDF[,1]),
                           Med=as.numeric(RectDF[,2]),
                           Max=as.numeric(RectDF[,3]),
                           Low = as.numeric(RectDF[,5]),
                           High = as.numeric(RectDF[,5]) + 1,
                           Stage = RectDF[,4])






# Gname <- sample(rownames(BuettProc[[1]]$NodesExp), 1)
Gname <- rownames(BuettProc[[1]]$NodesExp)[1]

tData <- NULL

for(i in 1:length(BuettProc)){
  
  if(!any(rownames(BuettProc[[i]]$NodesExp) == Gname))
    next()
  
  tData <- rbind(tData, cbind(BuettProc[[i]]$NodesPT/sum(BuettInfoList[[i]]$FinalStruct$ProjPoints[[length(BuettInfoList[[i]]$FinalStruct$ProjPoints)]]$EdgeLength),
                              BuettProc[[i]]$NodesExp[Gname, ],
                              rep(i, length(BuettProc[[i]]$NodesPT))))
}

colnames(tData) <- c("PT", "Exp", "Rep")
tDF <- data.frame(tData)
tDF$Rep <- as.character(tDF$Rep)


Loc.RectDF.Final <- RectDF.Final
Loc.RectDF.Final$Low <- Loc.RectDF.Final$Low - 1
Loc.RectDF.Final$High <- Loc.RectDF.Final$High - 1

MinExp <- min(tDF$Exp)
MaxExp <- max(tDF$Exp)
RangeExp <- MaxExp - MinExp
StepCat <- RangeExp/max(Loc.RectDF.Final$High)

Loc.RectDF.Final$Low <- Loc.RectDF.Final$Low*StepCat + MinExp
Loc.RectDF.Final$High <- Loc.RectDF.Final$High*StepCat + MinExp

p <- ggplot(data = tDF, mapping = aes(x = PT, y = Exp)) +
  geom_rect(data = Loc.RectDF.Final, mapping = aes(xmin = Min, xmax = Max, ymin = Low, ymax = High, fill = Stage),
            inherit.aes = FALSE, alpha = .3) + geom_smooth(span=.25) +
  geom_point(mapping = aes(color = Rep)) + labs(title = Gname)

print(p)



















RankVar <- apply(apply(ReordCells, 2, rank, ties.method = "min"), 1, var)
RankVar[which(RankVar > mean(RankVar))]

KM <- kmeans(RankVar, centers = range(RankVar))
if(KM$centers[1] < KM$centers[2]){
  ToUse <- which(KM$cluster == 1)
} else {
  ToUse <- which(KM$cluster == 2)
}













AllGenes <- lapply(BuettInfoList, function(x){
  colnames(x$FinalStruct$FiltExp)
})

AllGenesNales <- unique(unlist(AllGenes, use.names = FALSE))

GeneOrder <- sapply(AllGenes, function(x) {
  match(AllGenesNales, x)
})


AllGenesNales[which(rowSums(is.na(GeneOrder)) == 0)]


# AllFiltGenes <- rownames(GeneExprMat)
# 
# lapply(BuettInfoList, function(x){
#   AllFiltGenes <<- intersect(AllFiltGenes, x$Genes)
# })



FilteredFinalList <- list()
BuettProcList <- list()

for(i in 1:10){
  FilteredFinal <- ProjectAndCompute(DataSet = GeneExprMat,
                                     GeneSet = AllGenesNales[which(rowSums(is.na(GeneOrder)) == 0)],
                                     Categories = Categories, nNodes = 40,
                                     VarThr = .99, GraphType = 'Circle', PlanVarLimit = .85,
                                     PlanVarLimitIC = .9, ForceLasso = FALSE, InitStructNodes = 20,
                                     MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 20, DipPVThr = 1e-4,
                                     PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE, PCAFilter = TRUE,
                                     PlotDebug = FALSE)
  
  FilteredFinalList[[i]] <- FilteredFinal
  
  ProjectOnPrincipalGraph(Nodes = FilteredFinal$PrinGraph$Nodes, Edges = FilteredFinal$PrinGraph$Edges, Points = FilteredFinal$Data,
                          UsedPoints = NULL, Categories = FilteredFinal$Categories, Title=paste("(Filtered)"),
                          PCACenter = TRUE, ShowFitted = FALSE)
  
  BuettProc <- PlotOnStages(Structure = "Circle", TaxonList = FilteredFinal$TaxonList[[length(FilteredFinal$TaxonList)]],
                            Categories = FilteredFinal$Categories, nGenes = 2,
                            PrinGraph = FilteredFinal$PrinGraph,
                            Net = FilteredFinal$Net[[length(FilteredFinal$Net)]],
                            SelThr = .35, ComputeOverlaps = TRUE, ExpData = FilteredFinal$FiltExp,
                            RotatioMatrix = FilteredFinal$PCAData$rotation[,1:FilteredFinal$nDims],
                            PointProjections = FilteredFinal$ProjPoints[[length(FilteredFinal$ProjPoints)]])
  
  BuettProcList[[i]] <- BuettProc
}


















FilteredFinal <- ProjectAndCompute(DataSet = GeneExprMat,
                                   GeneSet = AllGenesNales[which(rowSums(is.na(GeneOrder)) == 0)],
                                   Categories = Categories, nNodes = 40,
                                   VarThr = .99, GraphType = 'Circle', PlanVarLimit = .85,
                                   PlanVarLimitIC = .9, ForceLasso = FALSE, InitStructNodes = 20,
                                   MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 20, DipPVThr = 1e-4,
                                   PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE, PCAFilter = TRUE,
                                   PlotDebug = FALSE)

ProjectOnPrincipalGraph(Nodes = FilteredFinal$PrinGraph$Nodes, Edges = FilteredFinal$PrinGraph$Edges, Points = FilteredFinal$Data,
                        UsedPoints = NULL, Categories = FilteredFinal$Categories, Title=paste("(Filtered)"),
                        PCACenter = TRUE, ShowFitted = FALSE)

BuettProc <- PlotOnStages(Structure = "Circle", TaxonList = FilteredFinal$TaxonList[[length(FilteredFinal$TaxonList)]],
                          Categories = FilteredFinal$Categories, nGenes = 2,
                          PrinGraph = FilteredFinal$PrinGraph,
                          Net = FilteredFinal$Net[[length(FilteredFinal$Net)]],
                          SelThr = .35, ComputeOverlaps = TRUE, ExpData = FilteredFinal$FiltExp,
                          RotatioMatrix = FilteredFinal$PCAData$rotation[,1:FilteredFinal$nDims],
                          PointProjections = FilteredFinal$ProjPoints[[length(FilteredFinal$ProjPoints)]])
















# Exapand genes

# WORK HERE!!!!



Topology = 'Circle'
DistillThr = .6
IgnoreTail = TRUE
Log = TRUE
StartQuant = .25
Title = "Buettner et al"
PlanVarLimit = .85
PlanVarLimitIC = .9
KeepOriginal = TRUE
PCACenter = FALSE
PlotDebug = FALSE
Mode = "VarPC"
ExtMode = 3
MaxRounds = 15




ProjectAndExpand <- function(GeneExprMat, StartSet, Categories, Topology = 'Circle', DistillThr = .6,
                             IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al",
                             PlanVarLimit = .85, PlanVarLimitIC = .9, KeepOriginal = TRUE,
                             PCACenter = FALSE, PlotDebug = FALSE, Mode = "VarPC", ExtMode = 3,
                             MaxRounds = 15, StopCrit = .95, ExpQuant = .01) {
  
  
  # Produce the initial analysis
  
  Info.Exp <- PlotsAndDistill(GeneExprMat = GeneExprMat, StartSet = StartSet,
                                         Categories = Categories, Topology = Topology, DistillThr = DistillThr,
                                         IgnoreTail = IgnoreTail, Log = Log, StartQuant = StartQuant, Title = Title, PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC,
                                         PCACenter = PCACenter, PlotDebug = PlotDebug, Mode = Mode, ExtMode = ExtMode)
  
  # Produce the initial processed data
  
  Proc.Exp <- PlotOnStages(Structure = Structure, TaxonList = Info.Exp$FinalStruct$TaxonList[[length(Info.Exp$FinalStruct$TaxonList)]],
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
  
  tictoc::tic()
  MedianVar <- apply(log10(GeneExprMat[AllMean >= min(SelMean), colnames(Proc.Exp$CellExp)] + 1), 1, function(x){
                                         AGG <- aggregate(x, by=list(TaxVect), median)
                                         AGGVect <- AGG[,2]
                                         names(AGGVect) <- paste(AGG[,1])
                                         median((aggregate(x, by=list(TaxVect), IQR))[,2]/AGGVect, na.rm=TRUE)
                                       })
  tictoc::toc()
  
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
  
  Round.Count <- 0
  
  while(!DONE){
    
    # Repeating the analysis untill convergence (or max iterations)
    
    Info.Exp.StageI <- PlotsAndDistill(GeneExprMat = GeneExprMat, StartSet = GeneStageList[[length(GeneStageList)]],
                                       Categories = Categories, Topology = Topology, DistillThr = DistillThr,
                                       IgnoreTail = IgnoreTail, Log = Log, StartQuant = StartQuant, Title = paste(Title, "STEP", Round.Count), PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC,
                                       PCACenter = PCACenter, PlotDebug = PlotDebug, Mode = Mode, ExtMode = ExtMode)
    
    Info.Exp.List[[length(Info.Exp.List)+1]] <- Info.Exp.StageI
    
    
    Proc.Exp.StageI <- PlotOnStages(Structure = "Circle", TaxonList = Info.Exp.StageI$FinalStruct$TaxonList[[length(Info.Exp.StageI$FinalStruct$TaxonList)]],
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
    
    Round.Count <- Round.Count + 1
    
    print(paste("Round", Round.Count, "Done"))
    
    if(Round.Count >= MaxRounds){
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
    
    tictoc::tic()
    MedianVar <- apply(log10(GeneExprMat[AllMean >= min(SelMean), colnames(Proc.Exp.StageI$CellExp)] + 1), 1, function(x){
      AGG <- aggregate(x, by=list(TaxVect), median)
      AGGVect <- AGG[,2]
      names(AGGVect) <- paste(AGG[,1])
      median((aggregate(x, by=list(TaxVect), IQR))[,2]/AGGVect, na.rm=TRUE)
    })
    tictoc::toc()
    
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
    
    if(GeneIntRatPrev > StopCrit & GeneIntRatCurr > StopCrit){
      print("Converged!")
      DONE <- TRUE
      break()
    }
    
  }
  
  return(list(StartInfo = Info.Exp, StartProc = Proc.Exp,
              ListInfo = Info.Exp.List, ListProc = Proc.Exp.List,
              GeneList = GeneStageList))
}






GeneExprMat <- Murine_ESC_All
Categories <- Murine_ESC_Stages





ExpDataset <- ProjectAndExpand(GeneExprMat = Murine_ESC_All, StartSet = MouseGenes_GOCellCycle, Categories = Murine_ESC_Stages,
                 Topology = 'Circle', DistillThr = .6, IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al",
                 PlanVarLimit = .85, PlanVarLimitIC = .9, KeepOriginal = TRUE, PCACenter = FALSE, PlotDebug = FALSE,
                 Mode = "VarPC", ExtMode = 3)
  





plot(unlist(lapply(ExpDataset$GeneList, length)), type='l')




length(intersect(ExpDataset$GeneList[[7]], ExpDataset$GeneList[[1]]))/length(ExpDataset$GeneList[[1]])

















BuettInfo.Exp.75 <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .75))), MouseGenes_GOCellCycle)),
                                          Categories = Categories, Topology = 'Circle', DistillThr = .5,
                                          IgnoreTail = TRUE, Log = TRUE, StartQuant = .75, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9,
                                          PCACenter = FALSE, PlotDebug = FALSE)

BuettInfo.Exp.25 <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .25))), MouseGenes_GOCellCycle)),
                                          Categories = Categories, Topology = 'Circle', DistillThr = .5,
                                          IgnoreTail = TRUE, Log = TRUE, StartQuant = .25, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9,
                                          PCACenter = FALSE, PlotDebug = FALSE)












BuettProc.Exp.25 <- PlotOnStages(Structure = "Circle", TaxonList = BuettInfo.Exp.25$FinalStruct$TaxonList[[length(BuettInfo.Exp.25$FinalStruct$TaxonList)]],
                                Categories = BuettInfo.Exp.25$FinalStruct$Categories, nGenes = 2,
                                PrinGraph = BuettInfo.Exp.25$FinalStruct$PrinGraph,
                                Net = BuettInfo.Exp.25$FinalStruct$Net[[length(BuettInfo.Exp.25$FinalStruct$Net)]],
                                SelThr = .35, ComputeOverlaps = TRUE, ExpData = BuettInfo.Exp.25$FinalStruct$FiltExp,
                                RotatioMatrix = BuettInfo.Exp.25$FinalStruct$PCAData$rotation[,1:BuettInfo.Exp.25$FinalStruct$nDims],
                                PointProjections = BuettInfo.Exp.25$FinalStruct$ProjPoints[[length(BuettInfo.Exp.25$FinalStruct$ProjPoints)]])


BuettProc.Exp.75 <- PlotOnStages(Structure = "Circle", TaxonList = BuettInfo.Exp.75$FinalStruct$TaxonList[[length(BuettInfo.Exp.75$FinalStruct$TaxonList)]],
                                 Categories = BuettInfo.Exp.75$FinalStruct$Categories, nGenes = 2,
                                 PrinGraph = BuettInfo.Exp.75$FinalStruct$PrinGraph,
                                 Net = BuettInfo.Exp.75$FinalStruct$Net[[length(BuettInfo.Exp.75$FinalStruct$Net)]],
                                 SelThr = .35, ComputeOverlaps = TRUE, ExpData = BuettInfo.Exp.75$FinalStruct$FiltExp,
                                 RotatioMatrix = BuettInfo.Exp.75$FinalStruct$PCAData$rotation[,1:BuettInfo.Exp.75$FinalStruct$nDims],
                                 PointProjections = BuettInfo.Exp.75$FinalStruct$ProjPoints[[length(BuettInfo.Exp.75$FinalStruct$ProjPoints)]])







barplot(rbind(c(length(BuettInfo.Exp$Genes), length(BuettInfo.Exp.25$Genes), length(BuettInfo.Exp.5$Genes), length(BuettInfo.Exp.75$Genes)),
              c(length(MouseGenes_GOCellCycle), length(unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .25))), MouseGenes_GOCellCycle))),
                length(unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .5))), MouseGenes_GOCellCycle))),
                length(unique(c(names(which(MedianVar < quantile(MedianVar[PrevDistilled], .75))), MouseGenes_GOCellCycle)))))
, beside = TRUE)













sum(BuettInfo$Genes %in% BuettInfo.Exp$Genes)/length(BuettInfo$Genes)
length(BuettInfo.Exp$Genes)/length(BuettInfo$Genes)




CellToPlot <- intersect(names(BuettProc.Exp$CellsPT), names(BuettProc$CellsPT))


plot(rank(BuettProc.Exp$CellsPT[CellToPlot], ties.method = "min"), rank(BuettProc$CellsPT[CellToPlot], ties.method = "min"))


cor(BuettProc.Exp$CellsPT[CellToPlot], BuettProc$CellsPT[CellToPlot])


cor(rank(BuettProc.Exp$CellsPT[CellToPlot], ties.method = "min"), rank(BuettProc$CellsPT[CellToPlot], ties.method = "min"))


GenesByStage <- lapply(BuettProc$StageOnNodes, function(x){
  BuettProc$NodesExp[,x]
})

StageOnNodes <- rep(NA, ncol(BuettProc$NodesExp))
lapply(as.list(1:length(BuettProc$StageOnNodes)),function(i){
  StageOnNodes[BuettProc$StageOnNodes[[i]]] <<- names(BuettProc$StageOnNodes)[i]
})
StageOnNodes <- factor(StageOnNodes, levels = names(BuettProc$StageOnNodes))

# boxplot(BuettProc$NodesExp[1,] ~ StageOnNodes)

AllPV <- apply(BuettProc$NodesExp, 1, function(x){
  AOV <- aov(x ~ StageOnNodes)
  summary(AOV)[[1]][1,5]
})


TLGenes <- apply(BuettProc$NodesExp[which(AllPV < 1e-3),], 1, function(x){
  AGG <- aggregate(x, by = list(StageOnNodes), median)
  TS <- AGG[which.max(AGG[,2]),1]
  LS <- AGG[which.min(AGG[,2]),1]

  c(as.character(TS), wilcox.test(x ~ StageOnNodes == TS, alternative = "less")$p.value,
    as.character(LS), wilcox.test(x ~ StageOnNodes == LS, alternative = "greater")$p.value)

})

TLGenes <- t(TLGenes)

table(TLGenes[TLGenes[,4] < 1e-3, 3])
table(TLGenes[TLGenes[,2] < 1e-3, 1])

# 
# 
# 
# AllPV <- apply(BuettProc$NodesExp, 1, function(x){
#   AOV <- aov(x ~ StageOnNodes)
#   PV <- summary(AOV)[[1]][1,5]
#   THT <- TukeyHSD(AOV)
#   SigDiff <- which(THT$StageOnNodes[,"p adj"] < 1e-3)
#   
#   RestTHT <- THT$StageOnNodes[SigDiff,]
#   SplGr <- strsplit(rownames(RestTHT), "-")
#   Gr1 <- unlist(lapply(SplGr, "[[", 1))
#   Gr2 <- unlist(lapply(SplGr, "[[", 2))
#   DifSgn <- sign(RestTHT[, "diff"])
# })
# 
# SelGenes <- which(AllPV < 1e-3)
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
