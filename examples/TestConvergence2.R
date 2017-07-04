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

FixGeneConv <- FALSE

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

if(!FixGeneConv){
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
  
} else {
  AllFreman_Mouse <- c(Freeman_CC1, Freeman_CC2, Freeman_CC9,
                       Freeman_G1S_CC4, Freeman_G1S_CC4A, Freeman_G2_CC6A,
                       Freeman_G2M_CC6, Freeman_M_CC6B, Freeman_S_CC4B)
  
  AllFreman_Mouse_CCP <- c(Freeman_G1S_CC4, Freeman_G1S_CC4A, Freeman_G2_CC6A,
                           Freeman_G2M_CC6, Freeman_M_CC6B, Freeman_S_CC4B)
  
  AllFreman_Mouse <- AllFreman_Mouse[!is.na(AllFreman_Mouse)]
  AllFreman_Mouse_CCP <- AllFreman_Mouse_CCP[!is.na(AllFreman_Mouse_CCP)]
  
  AllFreman_Mouse <- gsub("(\\b[a-z]{1})", "\\U\\1", tolower(AllFreman_Mouse), perl=TRUE)
  AllFreman_Mouse_CCP <- gsub("(\\b[a-z]{1})", "\\U\\1", tolower(AllFreman_Mouse_CCP), perl=TRUE)
  
  tGMT <- rRoma::SelectFromMSIGdb("GO_CELL_CYCLE")
  MouseGenes_GOCellCycle <- tGMT[[which(unlist(lapply(tGMT, "[[", "Name"))=="GO_CELL_CYCLE")]]$Genes
  
  MouseGenes_GOCellCycle <- gsub("(\\b[a-z]{1})", "\\U\\1", tolower(MouseGenes_GOCellCycle), perl=TRUE)
}









# MOUSE


# Define common functions ----------------------------------------------------------------


PlotsAndDistill.Mouse <- function(GeneExprMat, StartSet, Categories, Topology = "Circle", IgnoreTail = FALSE, PlanVarLimit = .9, PlanVarLimitIC = 92,
                                  DistillThr = .7, Log = TRUE, StartQuant = .5, Title = '', MinProlCells = 20) {
  
  PrGraph.Initial <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = StartSet, OutThr = 3, nNodes = 40,
                                       VarThr = .99, Categories = Categories, GraphType = Topology, PlanVarLimit = PlanVarLimit,
                                       PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, InitStructNodes = 20,
                                       MinBranDiff = 2, Log = Log, Filter = TRUE, MinProlCells = MinProlCells, DipPVThr = 1e-4,
                                       PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
  
  DistilledGenes <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = "VarPC", DistillThr = DistillThr, FastReduce = FALSE,
                                    ExtMode = 2, StopCrit = .95, OutThr = 3, nNodes = 40, VarThr = .99, Categories = Categories,
                                    GraphType = Topology, PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE,
                                    InitStructNodes = 20, MinBranDiff = 2, Log = Log, Filter = TRUE, MinProlCells = MinProlCells, DipPVThr = 1e-4,
                                    PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE,
                                    IgnoreTail = IgnoreTail, StartQuant = StartQuant)
  
  PrGraph.Final <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = DistilledGenes, OutThr = 3, nNodes = 40, 
                                     VarThr = .99, Categories = Categories, GraphType = Topology, PlanVarLimit = PlanVarLimit,
                                     PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, InitStructNodes = MinProlCells,
                                     MinBranDiff = 2, Log = Log, Filter = TRUE, MinProlCells = 20, DipPVThr = 1e-4,
                                     PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
  
  
  
  
  ProjectOnPrincipalGraph(Nodes = PrGraph.Initial$PrinGraph$Nodes, Edges = PrGraph.Initial$PrinGraph$Edges, Points = PrGraph.Initial$Data,
                          UsedPoints = NULL, Categories = PrGraph.Initial$Categories, Title=paste(Title, "(GO Cell Cycle)"),
                          PCACenter = TRUE, ShowFitted = FALSE)
  
  ProjectOnPrincipalGraph(Nodes = PrGraph.Final$PrinGraph$Nodes, Edges = PrGraph.Final$PrinGraph$Edges, Points = PrGraph.Final$Data,
                          UsedPoints = NULL, Categories = PrGraph.Final$Categories, Title=paste(Title, "(Filtered)"),
                          PCACenter = TRUE, ShowFitted = FALSE)
  
  if(!FixGeneConv){
    AllWit_Mouse <- gsub("(\\b[a-z]{1})", "\\U\\1", tolower(unlist(StageAssociation_Whit_Ext[-c(1,2)], use.names = FALSE)), perl=TRUE)
  } else {
    AllWit_Mouse <- ConvertNames("human", "mouse", unlist(StageAssociation_Whit_Ext[-c(1,2)], use.names = FALSE)) 
  }
  
  
  
  barplot(c(length(StartSet), length(DistilledGenes)), ylab = "Number of genes", names.arg = c("Base GO", "Filtered"))
  
  CIA <- prop.test(length(intersect(AllWit_Mouse, DistilledGenes)), length(DistilledGenes), conf.level = .95)$conf.int[1:2]*100
  CIB <- prop.test(length(intersect(AllWit_Mouse, StartSet)), length(StartSet), conf.level = .95)$conf.int[1:2]*100
  
  PT1 <- prop.test(c(length(intersect(AllWit_Mouse, DistilledGenes)), length(intersect(AllWit_Mouse, StartSet))),
                   c(length(DistilledGenes), length(StartSet)))
  
  B <- barplot(100*c(length(intersect(AllWit_Mouse, DistilledGenes)),
                     length(intersect(AllWit_Mouse, StartSet)))/c(length(DistilledGenes), length(StartSet)),
               names.arg = c("Filtered", "Base GO"), ylab = "Percentage of genes identified by Whithfield et al.", ylim = c(0, 25))
  
  arrows(x0 = B[1], x1 = B[1], y0 = CIA[1], y1 = CIA[2], angle = 90, length = .5, lwd = 2, code = 3)
  arrows(x0 = B[2], x1 = B[2], y0 = CIB[1], y1 = CIB[2], angle = 90, length = .5, lwd = 2, code = 3)
  
  
  CIA <- prop.test(length(intersect(AllFreman_Mouse, DistilledGenes)), length(DistilledGenes), conf.level = .95)$conf.int[1:2]*100
  CIB <- prop.test(length(intersect(AllFreman_Mouse, StartSet)), length(StartSet), conf.level = .95)$conf.int[1:2]*100
  
  PT2 <- prop.test(c(length(intersect(AllFreman_Mouse, DistilledGenes)), length(intersect(AllFreman_Mouse, StartSet))),
                   c(length(DistilledGenes), length(StartSet)))
  
  B <- barplot(100*c(length(intersect(AllFreman_Mouse, DistilledGenes)), 
                     length(intersect(AllFreman_Mouse, StartSet)))/c(length(DistilledGenes), length(StartSet)),
               names.arg = c("Filtered", "Base GO"), ylab = "Percentage of genes identified by Freeman et al.", ylim = c(0, 30))
  
  arrows(x0 = B[1], x1 = B[1], y0 = CIA[1], y1 = CIA[2], angle = 90, length = .5, lwd = 2, code = 3)
  arrows(x0 = B[2], x1 = B[2], y0 = CIB[1], y1 = CIB[2], angle = 90, length = .5, lwd = 2, code = 3)
  
  
  return(list(Genes = DistilledGenes, WitPT = PT1, FreePT = PT2, FinalStruct = PrGraph.Final))
  
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

Murine_ESC_Stages <- factor(Murine_ESC_Stages, levels = c("G1", "S", "G2M"))

GeneExprMat <- Murine_ESC_All
StartSet <- MouseGenes_GOCellCycle
Categories <- Murine_ESC_Stages

# library(BPSC)
# 
# 
# mat.res=estimateBPMatrix(data.matrix(GeneExprMat),para.num=4,fout=NULL,estIntPar=FALSE,useParallel=FALSE)
# 
# plot(apply(GeneExprMat, 1, mean), apply(GeneExprMat, 1, var), log="xy")
# 
# 
# 
# hist(as.integer(GeneExprMat[50,]))



# 
# # Filter data using scater
# 
# table(Categories, unlist(lapply(strsplit(colnames(GeneExprMat), "_"), "[[", 1)))
# 
# InfoData <- data.frame(Cell=colnames(GeneExprMat), Stage=Murine_ESC_Stages)
# rownames(InfoData) <- colnames(GeneExprMat)
# 
# library(scater)
# 
# pd <- new("AnnotatedDataFrame", InfoData)
# SceSet <- newSCESet(countData = GeneExprMat, phenoData = pd)
# 
# keep_features <- rowSums(exprs(SceSet) > 0) > 0
# SceSet <- SceSet[keep_features, ]
# 
# SceSet <- calculateQCMetrics(SceSet, feature_controls = 1:40)
# # scater_gui(SceSet)
# plotQC(SceSet, type = "exprs")







BuettInfo <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .6,
                                   IgnoreTail = TRUE, Log = TRUE, StartQuant = .55, Title = "Buettner et al", PlanVarLimit = .85, PlanVarLimitIC = .9)


PlotOnStages(Structure = "Circle", TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]],
             Categories = BuettInfo$FinalStruct$Categories, nGenes = 50,
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
# nGenes = 2
# PrinGraph = BuettInfo$FinalStruct$PrinGraph
# Net = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]]
# SelThr = .35
# ComputeOverlaps = TRUE
# ExpData = BuettInfo$FinalStruct$FiltExp
# RotatioMatrix = BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims]
# PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]]









# 
# 
# GeneExprMat = GeneExprMat
# StartSet = StartSet
# Categories = Categories
# Topology = 'Circle'
# DistillThr = .6
# IgnoreTail = TRUE
# Log = TRUE
# StartQuant = .55
# Title = "Buettner et al"
# PlanVarLimit = .85
# PlanVarLimitIC = .9



# 
# 
# 
# BuettInfo$WitPT
# BuettInfo$FreePT
# 
# TaxProj <- BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]]
# 
# TaxVect <- rep(NA, max(unlist(TaxProj), na.rm = TRUE))
# for(i in 1:length(TaxProj)){
#   TaxVect[TaxProj[[i]]] <- i 
# }
# 
# TB <- table(BuettInfo$FinalStruct$Categories, TaxVect)
# colnames(TB)
# 
# AllPaths <- GetLongestPath(Net = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
#                            Structure = "Circle", Circular = TRUE)
# 
# SummInfo <- NULL
# 
# for(i in 1:nrow(AllPaths$VertNumb)){
#   SelPath <- AllPaths$VertNumb[i,]
#   SelPath <- SelPath[SelPath %in% unique(TaxVect)]
#   SelPath <- SelPath[-length(SelPath)]
#   
#   Reordered <- TaxVect
#   for(j in 1:length(SelPath)){
#     Reordered[TaxVect == SelPath[j]] <- j
#   }
#   
#   AGG <- aggregate(Reordered, by = list(BuettInfo$FinalStruct$Categories), median)
#   AGG2 <- aggregate(Reordered, by = list(BuettInfo$FinalStruct$Categories), min)
#   AGG3 <- aggregate(Reordered, by = list(BuettInfo$FinalStruct$Categories), max)
#   
#   
#   if(isSorted(AGG[,2])){
#     SummInfo <- rbind(SummInfo,
#                       c(i, 1, summary(aov(Reordered ~ BuettInfo$FinalStruct$Categories))[[1]][1,"Pr(>F)"],
#                         AGG2[1,2], AGG[1,2])
#     )
#     
#     boxplot(Reordered ~ BuettInfo$FinalStruct$Categories, main = i)
#     
#   }
#   
#   
#   if(isSorted(rev(AGG[,2]))){
#     SummInfo <- rbind(SummInfo,
#                       c(i, 2, summary(aov(Reordered ~ BuettInfo$FinalStruct$Categories))[[1]][1,"Pr(>F)"],
#                         nrow(BuettInfo$FinalStruct$PrinGraph$Nodes) - AGG3[1,2] + 1,
#                         nrow(BuettInfo$FinalStruct$PrinGraph$Nodes) - AGG[1,2] + 1)
#     )
#     boxplot(Reordered ~ BuettInfo$FinalStruct$Categories, main = paste(i, "rev"))
#   }
# }
# 
# 
# 
# 
# Selected <- SummInfo[SummInfo[, 3] == min(SummInfo[, 3]), ]
# 
# if(length(Selected)>5){
#   Selected <- Selected[which.min(Selected[,4]),]
# }
# 
# SelPath <- AllPaths$VertNumb[Selected[1],]
# SelPathBuett <- SelPath
# 
# SelPath <- SelPath[SelPath %in% unique(TaxVect)]
# SelPath <- SelPath[-length(SelPath)]
# 
# if(Selected[2] == 2){
#   SelPath <- rev(SelPath)
#   SelPathBuett <- rev(SelPathBuett)
# }
# 
# Reordered <- TaxVect
# for(j in 1:length(SelPath)){
#   Reordered[TaxVect == SelPath[j]] <- j
# }
# 
# boxplot(Reordered ~ BuettInfo$FinalStruct$Categories)
# 
# ReorderedBuett <- Reordered




# Sasagawa et al ----------------------------------------------------------------


BaseDir <- "~/Google Drive/Datasets/Sasagawa et al - Murine Stem Cells (EB5 cell line)/"

FilesToRead <- list.files(path = paste(BaseDir, "ES/", sep = ''), full.names = TRUE, pattern = "txt.gz")
LoadedData <- list()

FilesToReadShort <- list.files(path = paste(BaseDir, "ES/", sep = ''), full.names = FALSE, pattern = "txt.gz")

ID <- unlist(lapply(strsplit(FilesToRead, "_"), "[[", 2))
StageVect <- rep("S_G1", length(ID))
StageVect[ID == "ESS"] <- "S_S"
StageVect[ID == "ESM"] <- "S_G2/M"

StageVect <- factor(StageVect, levels = c("S_G1", "S_S", "S_G2/M"))

for(i in 1:length(FilesToRead)){
  Data <- read_delim(FilesToRead[i], "\t", escape_double = FALSE, trim_ws = TRUE)
  LoadedData[[i]] <- list(Name = FilesToRead[i], Data = Data, Stage = StageVect[i])
}

AllGenes <- NULL

for(i in 1:length(LoadedData)){
  AllGenes <- unique(c(AllGenes, LoadedData[[i]]$Data$gene.symbol))
}

AllGenesEns <- rep(NA, length(AllGenes))
names(AllGenesEns) <- AllGenes

for(i in 1:length(LoadedData)){
  AllGenesEns[LoadedData[[i]]$Data$gene.symbol] <- LoadedData[[i]]$Data$id
}

ValMat <- matrix(rep(0, length(AllGenes)*length(FilesToRead)), ncol = length(AllGenes))

colnames(ValMat) <- AllGenes
rownames(ValMat) <- FilesToRead

for(i in 1:length(LoadedData)){
  ValMat[i,LoadedData[[i]]$Data$gene.symbol] <- LoadedData[[i]]$Data$fpkm
}


GeneExprMat <- t(data.matrix(ValMat))
StartSet <- MouseGenes_GOCellCycle
Categories <- factor(StageVect)


SasaInfo <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .6,
                                  IgnoreTail = TRUE, Log = TRUE, StartQuant = .55, Title = "Sasagawa et al", PlanVarLimit = .85, PlanVarLimitIC = .9)


SasaInfo$WitPT
SasaInfo$FreePT





TaxProj <- SasaInfo$FinalStruct$TaxonList[[length(SasaInfo$FinalStruct$TaxonList)]]

TaxVect <- rep(NA, max(unlist(TaxProj), na.rm = TRUE))
for(i in 1:length(TaxProj)){
  TaxVect[TaxProj[[i]]] <- i 
}

TB <- table(SasaInfo$FinalStruct$Categories, TaxVect)
colnames(TB)

AllPaths <- GetLongestPath(Net = SasaInfo$FinalStruct$Net[[length(SasaInfo$FinalStruct$Net)]],
                           Structure = "Circle", Circular = TRUE)


SummInfo <- NULL

for(i in 1:nrow(AllPaths$VertNumb)){
  SelPath <- AllPaths$VertNumb[i,]
  SelPath <- SelPath[SelPath %in% unique(TaxVect)]
  SelPath <- SelPath[-length(SelPath)]
  
  Reordered <- TaxVect
  for(j in 1:length(SelPath)){
    Reordered[TaxVect == SelPath[j]] <- j
  }
  
  # boxplot(Reordered ~ SasaInfo$FinalStruct$Categories)
  
  AGG <- aggregate(Reordered, by = list(SasaInfo$FinalStruct$Categories), median)
  
  if(isSorted(AGG[,2])){
    SummInfo <- rbind(SummInfo,
                      c(i, 1, summary(aov(Reordered ~ SasaInfo$FinalStruct$Categories))[[1]][1,"Pr(>F)"])
    )
  }
  
  
  if(isSorted(rev(AGG[,2]))){
    SummInfo <- rbind(SummInfo,
                      c(i, 2, summary(aov(Reordered ~ SasaInfo$FinalStruct$Categories))[[1]][1,"Pr(>F)"])
    )
  }
}

Selected <- SummInfo[which.min(SummInfo[,3]), ]

SelPath <- AllPaths$VertNumb[Selected[1],]
SelPathSasa <- SelPath

SelPath <- SelPath[SelPath %in% unique(TaxVect)]
SelPath <- SelPath[-length(SelPath)]

if(Selected[2] == 2){
  SelPath <- rev(SelPath)
  SelPathSasa <- rev(SelPathSasa)
}

Reordered <- TaxVect
for(j in 1:length(SelPath)){
  Reordered[TaxVect == SelPath[j]] <- j
}

boxplot(Reordered ~ SasaInfo$FinalStruct$Categories)

ReorderedSasa <- Reordered
















# Buettner et al VS Sasagawa et al ----------------------------------------------------------------

pie(c(length(setdiff(BuettInfo$Genes, SasaInfo$Genes)),
      length(intersect(BuettInfo$Genes, SasaInfo$Genes)),
      length(setdiff(SasaInfo$Genes, BuettInfo$Genes))),
    labels = c("Buettner only", "Shared", "Sasagawa only"),
    main = "Gene filtering overlap")





Group1 <- lapply(as.list(1:100), sample, x = MouseGenes_GOCellCycle, size = length(BuettInfo$Genes))
Group2 <- lapply(as.list(1:100), sample, x = MouseGenes_GOCellCycle, size = length(SasaInfo$Genes))
IntersSamp <- unlist(lapply(mapply(intersect, x = Group1, y = Group2), length))
IntersReal <- length(intersect(BuettInfo$Genes, SasaInfo$Genes))

boxplot(IntersSamp, ylim = c(min(IntersSamp)-10, max(IntersReal)+10), at = 1, ylab = "Number of shared genes")
points(x = 1, y =IntersReal, cex = 3, col="red", pch=20)












# Buettner et al and Sasagawa et al peaks ----------------------------------------------------------------

# Structure <- "Circle"
# TaxonList <- BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]]
# Categories <- 
# Net <- BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]]
# RotatioMatrix <- BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims]
# PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]]
# ExpData <- BuettInfo$FinalStruct$FiltExp


# 
# TaxProj <- BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]]
# Empty <- which(unlist(lapply(lapply(TaxProj, is.na), any)))
# 
# TaxVect <- rep(NA, max(unlist(TaxProj), na.rm = TRUE))
# for(i in 1:length(TaxProj)){
#   TaxVect[TaxProj[[i]]] <- i 
# }
# 
# TB <- table(BuettInfo$FinalStruct$Categories, TaxVect)
# colnames(TB)
# 
# 
# 
# ExtendedTB <- TB
# 
# for(i in 1:length(Empty)){
#   ExtendedTB <- cbind(ExtendedTB, rep(0, length(levels(BuettInfo$FinalStruct$Categories))))
#   colnames(ExtendedTB)[ncol(ExtendedTB)] <- paste(Empty[i]) 
# }
# 
# 
# # AllPaths <- GetLongestPath(Net = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
# #                            Structure = "Circle", Circular = TRUE)
# # 
# # SummMediStat <- NULL
# # 
# # for(i in 1:nrow(AllPaths$VertNumb)){
# #   SelPath <- AllPaths$VertNumb[i,]
# #   SelPath <- SelPath[SelPath %in% unique(TaxVect)]
# #   
# #   ExtendedTB <- ExtendedTB[,SelPath]
# #   
# #   S1 <- unlist(mapply(rep, 1:ncol(ExtendedTB), ExtendedTB[1,]))
# #   S2 <- unlist(mapply(rep, 1:ncol(ExtendedTB), ExtendedTB[2,]))
# #   S3 <- unlist(mapply(rep, 1:ncol(ExtendedTB), ExtendedTB[3,]))
# #   
# #   SummMediStat <- rbind(SummMediStat, 
# #                         c(
# #                           kruskal.test(c(S1, S2, S3), factor(c(rep("S1", length(S1)), rep("S2", length(S2)), rep("S3", length(S3)))))$p.value,
# #                           median(S1),
# #                           median(S3)
# #                         )
# #   )
# # }
# 
# 
# # SelPath <- AllPaths$VertNumb[which.min(SummMediStat[,1]),]
# # 
# # if(!isSorted(SummMediStat[which.min(SummMediStat[,1]),2:3])){
# #   SelPath <- rev(SelPath)
# # }
# 
# # SelPath <- SelPath[SelPath %in% unique(TaxVect)]
# ExtendedTB <- ExtendedTB[,SelPathBuett]
# 
# barplot(t(t(ExtendedTB)/colSums(ExtendedTB)), col = c("red", "green", "blue"), beside = TRUE, las = 2)
# 
# 
# 
# SelPath <- SelPath[-length(SelPath)]
# SelPathSTG <- rep(NA, length(SelPath))
# 
# 
# PercMat <- t(t(ExtendedTB)/colSums(ExtendedTB))
# PercMat[is.na(PercMat)] <- 0
# BinPercMat <- (PercMat > .5)
# 
# if(!any(BinPercMat[,1])){
#   BinPercMat[,1] <- BinPercMat[,2]
# }
# 
# 
# 
# 
# 
# 
# for(j in 1:nrow(BinPercMat)){
#   if(BinPercMat[j,2] & BinPercMat[j,ncol(BinPercMat)]){
#     BinPercMat[j,1] <- TRUE
#   }
# }
# 
# for(i in 2:(ncol(BinPercMat)-1)){
#   for(j in 1:nrow(BinPercMat)){
#     if(BinPercMat[j,i-1] & BinPercMat[j,i+1]){
#       BinPercMat[j,i] <- TRUE
#     }
#   }
# }
# 
# for(j in 1:nrow(BinPercMat)){
#   if(BinPercMat[j,1] & BinPercMat[j,ncol(BinPercMat)-1]){
#     BinPercMat[j,ncol(BinPercMat)] <- TRUE
#   }
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
# for(j in 1:nrow(BinPercMat)){
#   if(!BinPercMat[j,2] & !BinPercMat[j,ncol(BinPercMat)]){
#     BinPercMat[j,1] <- FALSE
#   }
# }
# 
# for(i in 2:(ncol(BinPercMat)-1)){
#   for(j in 1:nrow(BinPercMat)){
#     if(!BinPercMat[j,i-1] & !BinPercMat[j,i+1]){
#       BinPercMat[j,i] <- FALSE
#     }
#   }
# }
# 
# for(j in 1:nrow(BinPercMat)){
#   if(!BinPercMat[j,1] & !BinPercMat[j,ncol(BinPercMat)-1]){
#     BinPercMat[j,ncol(BinPercMat)] <- FALSE
#   }
# }
# 
# 
# 
# 
# 
# G1S <- BinPercMat[1,] & BinPercMat[2,]
# 
# SG2 <- BinPercMat[2,] & BinPercMat[3,]
# 
# MG1 <- BinPercMat[1,] & BinPercMat[3,]
# 
# BinPercMatExt <- rbind(BinPercMat[1,] & !G1S & !MG1, G1S,
#                        BinPercMat[2,] & !G1S & !SG2, SG2,
#                        BinPercMat[3,] & !SG2 & !MG1, MG1)*1
# rownames(BinPercMatExt) <- c("G1", "G1S", "S", "SG2", "G2M", "MG1")
# 
# # while(all(BinPercMatExt[,1] == BinPercMatExt[,ncol(BinPercMatExt)])){
# #   BinPercMatExt <- cbind(BinPercMatExt[,1], BinPercMatExt[,1:(ncol(BinPercMatExt)-1)])
# # }
# 
# SelPathBuett <- SelPathBuett[-length(SelPathBuett)]
# BinPercMatExt <- BinPercMatExt[, -ncol(BinPercMatExt)]
# 
# while(any(BinPercMatExt[1:2,ncol(BinPercMatExt)]==1) & all(BinPercMatExt[3:6,ncol(BinPercMatExt)]!=1)){
#   BinPercMatExt <- cbind(BinPercMatExt[, ncol(BinPercMatExt)], BinPercMatExt[, -ncol(BinPercMatExt)])
#   SelPathBuett <- c(SelPathBuett[length(SelPathBuett)], SelPathBuett[-length(SelPathBuett)])
# }
# 
# while(any(BinPercMatExt[5:6,1]==1) & all(BinPercMatExt[1:4,1]!=1)){
#   BinPercMatExt <- cbind(BinPercMatExt[, -1], BinPercMatExt[, 1])
#   SelPathBuett <- c(SelPathBuett[-1], SelPathBuett[1])
# }
# 
# SelPathBuett <- c(SelPathBuett, SelPathBuett[1])
# BinPercMatExt <- cbind(BinPercMatExt, rep(0, 6))
# 
# 
# Idxs <- apply(BinPercMatExt==1, 1, which)
# 
# Bond <- lapply(Idxs[lapply(Idxs, length) > 1], range)
# 
# if(!isSorted(unlist(Bond))){
#   warning("Inconsitent data! Stages are not sequential!")
# }
# 
# 
# NodeOnGenes <- t(BuettInfo$FinalStruct$PrinGraph$Nodes %*% t(BuettInfo$FinalStruct$PCAData$rotation[,1:BuettInfo$FinalStruct$nDims]))
# 
# 
# OrderedPoints <- OrderOnPath(PrinGraph = BuettInfo$FinalStruct$PrinGraph, Path = SelPathBuett,
#                              PointProjections = BuettInfo$FinalStruct$ProjPoints[[length(BuettInfo$FinalStruct$ProjPoints)]])
# 
# 
# 
# # SelPath <- AllPaths$VertNumb[25,]
# # SelPath <- SelPath[-length(SelPath)]
# 
# 
# 
# 
# 
# 
# for(Idx in sample(1:nrow(NodeOnGenes), 25)){
#   
#   p <- ggplot(data = data.frame(x=cumsum(OrderedPoints$PathLen),
#                                 y=NodeOnGenes[Idx,as.numeric(SelPathBuett)]),
#               mapping = aes(x = x, y = y, color="PC")) + labs(x = "Pseudotime", y="Gene expression (log10 Pseudocount)", title = rownames(NodeOnGenes)[Idx]) +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   RangeDF <- data.frame(t(as.data.frame(lapply(Idxs, range))))
#   colnames(RangeDF) <- c("Min", "Max")
#   RangeDF[!is.finite(RangeDF$Min),] <- NA
#   RangeDF <- cbind(rownames(RangeDF), RangeDF)
#   colnames(RangeDF) <- c("Stage", "Min", "Max")
#   
#   RangeDF$Stage <- factor(as.character(RangeDF$Stage), levels = c("G1", "G1S", "S", "SG2", "G2M", "MG1"))
#   RangeDF$Min <- as.numeric(as.character(RangeDF$Min))
#   RangeDF$Max <- as.numeric(as.character(RangeDF$Max))
#   
#   RangeDF$Min <- cumsum(OrderedPoints$PathLen)[RangeDF$Min]
#   RangeDF$Max <- cumsum(OrderedPoints$PathLen)[RangeDF$Max]
#   
#   for(i in 2:nrow(RangeDF)){
#     RangeDF$Max[i-1] <- mean(c(RangeDF$Max[i-1], RangeDF$Min[i]), na.rm=TRUE)
#     RangeDF$Min[i] <- RangeDF$Max[i-1]
#   }
#   
#   p <- p + geom_rect(data = RangeDF, mapping = aes(fill=Stage, xmin=Min, xmax=Max),
#                      ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
#     geom_point(data = data.frame(x=OrderedPoints$PositionOnPath, y=BuettInfo$FinalStruct$FiltExp[,Idx]),
#                mapping = aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
#     geom_point() + geom_line() + scale_color_manual(values = c("blue", "black"))
#   
#   print(p)
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
# Idx <- 10
# 
# p <- ggplot(data = data.frame(x=1:(ncol(NodeOnGenes)+1),
#                          y=NodeOnGenes[Idx,as.numeric(SelPathBuett)]),
#        mapping = aes(x = x, y = y)) + labs(x = "Pseudotime", y="Gene expression (log10 Pseudocount)", title = rownames(NodeOnGenes)[Idx]) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# RangeDF <- data.frame(t(as.data.frame(lapply(Idxs, range))))
# colnames(RangeDF) <- c("Min", "Max")
# RangeDF[!is.finite(RangeDF$Min),] <- NA
# RangeDF <- cbind(rownames(RangeDF), RangeDF)
# colnames(RangeDF) <- c("Stage", "Min", "Max")
# 
# RangeDF$Stage <- factor(as.character(RangeDF$Stage), levels = c("G1", "G1S", "S", "SG2", "G2M", "MG1"))
# RangeDF$Min <- as.numeric(as.character(RangeDF$Min)) - .5
# RangeDF$Max <- as.numeric(as.character(RangeDF$Max)) + .5
# 
# 
# p <- p + geom_rect(data = RangeDF, mapping = aes(fill=Stage, xmin=Min, xmax=Max),
#               ymin = -Inf, ymax = Inf, inherit.aes = FALSE) + geom_point() + geom_line() 
# 
# print(p)




# WORK HEREEEEEEEEE




# 
# 
# TaxProj <- SasaInfo$FinalStruct$TaxonList[[length(SasaInfo$FinalStruct$TaxonList)]]
# 
# TaxVect <- rep(NA, max(unlist(TaxProj), na.rm = TRUE))
# for(i in 1:length(TaxProj)){
#   TaxVect[TaxProj[[i]]] <- i 
# }
# 
# 
# 
# 
# AllPaths <- GetLongestPath(Net = SasaInfo$FinalStruct$Net[[length(SasaInfo$FinalStruct$Net)]],
#                            Structure = "Circle", Circular = TRUE)
# SelPath <- AllPaths$VertNumb[6,]
# SelPath <- SelPath[SelPath %in% unique(TaxVect)]
# 
# TB <- table(SasaInfo$FinalStruct$Categories, as.character(TaxVect))
# TB <- TB[,rev(SelPath)]
# TB
# 
# barplot(t(t(TB)/colSums(TB)), col = c("red", "green", "blue"), beside = TRUE, las = 2)
# 
# SelPath <- SelPath[-length(SelPath)]
# SelPathSTG <- rep(NA, length(SelPath))
# 
# G1 <- SelPath[1:3]
# S <- SelPath[4:5]
# G2 <- SelPath[6:length(SelPath)]
# 
# SelPathSTG[1:3] <- "G1"
# SelPathSTG[4:5] <- "S"
# SelPathSTG[6:length(SelPath)] <- "G2"
# 
# 
# 
# 
# 
# STGVect <- TaxVect
# STGVect[TaxVect %in% G1] <- "G1"
# STGVect[TaxVect %in% S] <- "S"
# STGVect[TaxVect %in% G2] <- "G2"
# 
# 
# GeneDiff <- apply(BuettInfo$FinalStruct$FiltExp, 2, kruskal.test, g = factor(STGVect))
# GeneDiffPV <- unlist(lapply(GeneDiff, "[[", "p.value"))
# 
# aggregate(x = BuettInfo$FinalStruct$FiltExp[,1], by = list(STGVect), FUN = median)
# 
# AggList <- apply(X = BuettInfo$FinalStruct$FiltExp, MARGIN = 2, FUN = aggregate, by = list(STGVect), median)
# 
# MaxList <- unlist(lapply(AggList, function(x){which.max(x[, 2])}))
# 
# G1_Genes_Buett <- intersect(names(which(MaxList == 1)), names(GeneDiffPV[GeneDiffPV < 5e-2]))
# S_Genes_Buett <- intersect(names(which(MaxList == 4)), names(GeneDiffPV[GeneDiffPV < 5e-2]))
# G2M_Genes_Buett <- intersect(names(which(MaxList == 3)), names(GeneDiffPV[GeneDiffPV < 5e-2]))
# G1S_Genes_Buett <- intersect(names(which(MaxList == 2)), names(GeneDiffPV[GeneDiffPV < 5e-2]))
# 
# 
# StageAssociation_Buett_Inferred <- list(QVarCutOff = .75, Stages = c("G1", "S", "G2M"), S1_U = G1_Genes_Buett,
#                                         S2_U = G1S_Genes_Buett, S3_U = S_Genes_Buett, S4_U = G2M_Genes_Buett)
# 
# 
# 
# intersect(StageAssociation_Whit$S1_U, toupper(G1_Genes_Buett))
# intersect(StageAssociation_Whit$S2_U, toupper(S_Genes_Buett))
# intersect(StageAssociation_Whit$S3_U, toupper(G2M_Genes_Buett))
# intersect(StageAssociation_Whit$S4_U, toupper(G2M_Genes_Buett))
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






















GeneDiff <- apply(SasaInfo$FinalStruct$FiltExp, 2, kruskal.test, g = factor(STGVect))
GeneDiffPV <- unlist(lapply(GeneDiff, "[[", "p.value"))

aggregate(x = SasaInfo$FinalStruct$FiltExp[,1], by = list(STGVect), FUN = median)

AggList <- apply(X = SasaInfo$FinalStruct$FiltExp, MARGIN = 2, FUN = aggregate, by = list(STGVect), median)

MaxList <- unlist(lapply(AggList, function(x){which.max(x[, 2])}))

G1_Genes_Sasa <- intersect(names(which(MaxList == 1)), names(GeneDiffPV[GeneDiffPV < 5e-2]))
S_Genes_Sasa <- intersect(names(which(MaxList == 3)), names(GeneDiffPV[GeneDiffPV < 5e-2]))
G2M_Genes_Sasa <- intersect(names(which(MaxList == 2)), names(GeneDiffPV[GeneDiffPV < 5e-2]))


StageAssociation_Sasa_Inferred <- list(QVarCutOff = .75, Stages = c("G1", "S", "G2M"), S1_U = G1_Genes_Sasa,
                                       S2_U = S_Genes_Sasa, S3_U = G2M_Genes_Sasa)



intersect(StageAssociation_Whit$S1_U, toupper(G1_Genes_Sasa))
intersect(StageAssociation_Whit$S2_U, toupper(S_Genes_Sasa))
intersect(StageAssociation_Whit$S3_U, toupper(G2M_Genes_Sasa))
intersect(StageAssociation_Whit$S4_U, toupper(G2M_Genes_Sasa))




length(intersect(G1_Genes_Buett, G1_Genes_Sasa))/length(G1_Genes_Buett)
length(intersect(S_Genes_Buett, S_Genes_Sasa))/length(S_Genes_Buett)
length(intersect(G2M_Genes_Buett, G2M_Genes_Sasa))/length(G2M_Genes_Buett)




save(StageAssociation_Sasa_Inferred, StageAssociation_Buett_Inferred, BuettInfo, SasaInfo,
     file = "~/Google Drive/Datasets/Gene List/InferredGenes.RData")
















# Kowalczyk et al (C57BL6) ----------------------------------------------------------------

library(readxl)

BaseDir <- "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/"

Murine_HSC_Data <- read_excel(paste(BaseDir, "GSE59114_C57BL6_GEO_all.xlsx", sep = ''))
Murine_HSC_Samples <- read_excel(paste(BaseDir, "Table_S2.xlsx", sep = ''))

CCVect <- factor(Murine_HSC_Samples$`Estimated phase`, levels = c("G0", "G1(early)", "G1(late)", "S", "G2/M"))
names(CCVect) <- Murine_HSC_Samples$`cell ID`

ProgVect <- Murine_HSC_Samples$`Progression Rank`
names(ProgVect) <- Murine_HSC_Samples$`cell ID`

SelectedSamples <- Murine_HSC_Samples$`cell ID`[Murine_HSC_Samples$population %in% c("'old LT-HSC'", "'old LT-HSC(rep)'", "'old ST-HSC'", 
                                                                                     "'old ST-HSC(rep)'", "'young LT-HSC'", "'young ST-HSC'",
                                                                                     "'young MPP'", "'old MPP'")]
SelectedSamples <- intersect(SelectedSamples, colnames(Murine_HSC_Data))

Categories <- CCVect[paste(SelectedSamples)]

AllData <- data.matrix(Murine_HSC_Data[-1,paste(SelectedSamples)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories

KowalczykInfo.C57BL6.All <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .5,
                                                  IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al (All)",
                                                  PlanVarLimit = .85, PlanVarLimitIC = .9)

KowalczykInfo.C57BL6.All$WitPT
KowalczykInfo.C57BL6.All$FreePT



PlotOnStages(Structure = "Circle", TaxonList = KowalczykInfo.C57BL6.All$FinalStruct$TaxonList[[length(KowalczykInfo.C57BL6.All$FinalStruct$TaxonList)]],
             Categories = KowalczykInfo.C57BL6.All$FinalStruct$Categories, nGenes = 10,
             PrinGraph = KowalczykInfo.C57BL6.All$FinalStruct$PrinGraph,
             Net = KowalczykInfo.C57BL6.All$FinalStruct$Net[[length(KowalczykInfo.C57BL6.All$FinalStruct$Net)]],
             SelThr = .35, ComputeOverlaps = TRUE, ExpData = KowalczykInfo.C57BL6.All$FinalStruct$FiltExp,
             RotatioMatrix = KowalczykInfo.C57BL6.All$FinalStruct$PCAData$rotation[,1:KowalczykInfo.C57BL6.All$FinalStruct$nDims],
             PointProjections = KowalczykInfo.C57BL6.All$FinalStruct$ProjPoints[[length(KowalczykInfo.C57BL6.All$FinalStruct$ProjPoints)]])



ColCat <- c("red", "blue", "green")
names(ColCat) <- c("G1", "S", "G2M")


plotPieNet(Results =  BuettInfo$FinalStruct$IntGrahs[[length(BuettInfo$FinalStruct$IntGrahs)]],
           Data = BuettInfo$FinalStruct$Data, NodeSizeMult = 4,
           Categories = BuettInfo$FinalStruct$Categories, PlotNet = TRUE,
           Graph = BuettInfo$FinalStruct$Net[[length(BuettInfo$FinalStruct$Net)]],
           TaxonList = BuettInfo$FinalStruct$TaxonList[[length(BuettInfo$FinalStruct$TaxonList)]],
           LayOut = 'circle', Main = "Pincipal graph", ColCat = ColCat)

legend(x = "center", legend = names(ColCat), fill = ColCat)















KowalczykInfo.C57BL6.All.Tail <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Lasso', DistillThr = .5,
                                                  IgnoreTail = TRUE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al (All)",
                                                  PlanVarLimit = .85, PlanVarLimitIC = .9, MinProlCells = 100)



SelectedSamples <- Murine_HSC_Samples$`cell ID`[Murine_HSC_Samples$population %in% c("'old LT-HSC'", "'old LT-HSC(rep)'", "'young LT-HSC'")]
SelectedSamples <- intersect(SelectedSamples, colnames(Murine_HSC_Data))

Categories <- CCVect[paste(SelectedSamples)]

AllData <- data.matrix(Murine_HSC_Data[-1,paste(SelectedSamples)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories

KowalczykInfo.C57BL6.LTHSC <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .5,
                                                    IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al (LT-HSC)", PlanVarLimit = .85, PlanVarLimitIC = .9)

KowalczykInfo.C57BL6.LTHSC$WitPT
KowalczykInfo.C57BL6.LTHSC$FreePT





SelectedSamples <- Murine_HSC_Samples$`cell ID`[Murine_HSC_Samples$population %in% c("'old ST-HSC'", "'old ST-HSC(rep)'", "'young ST-HSC'")]
SelectedSamples <- intersect(SelectedSamples, colnames(Murine_HSC_Data))

Categories <- CCVect[paste(SelectedSamples)]

AllData <- data.matrix(Murine_HSC_Data[-1,paste(SelectedSamples)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories

KowalczykInfo.C57BL6.STHSC <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .6,
                                                    IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al (ST-HSC)", PlanVarLimit = .85, PlanVarLimitIC = .85)

KowalczykInfo.C57BL6.STHSC$WitPT
KowalczykInfo.C57BL6.STHSC$FreePT




SelectedSamples <- Murine_HSC_Samples$`cell ID`[Murine_HSC_Samples$population %in% c("'young MPP'", "'old MPP'")]
SelectedSamples <- intersect(SelectedSamples, colnames(Murine_HSC_Data))

Categories <- CCVect[paste(SelectedSamples)]

AllData <- data.matrix(Murine_HSC_Data[-1,paste(SelectedSamples)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories

KowalczykInfo.C57BL6.MPP <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Circle', DistillThr = .6,
                                                  IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al (MPP)", PlanVarLimit = .85, PlanVarLimitIC = .9)

KowalczykInfo.C57BL6.MPP$WitPT
KowalczykInfo.C57BL6.MPP$FreePT















# Kowalczyk et al (DBA) / Load Data ----------------------------------------------------------------

BaseDir <- "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/"

Murine_HSC_Data <- read_excel(paste(BaseDir, "GSE59114_DBA_GEO_all.xlsx", sep = ''))

BulkId <- grep("population", unlist(Murine_HSC_Data[1,]), ignore.case = TRUE)
Murine_HSC_Data <- Murine_HSC_Data[, -BulkId]

Young <- grep("young", unlist(Murine_HSC_Data[1,]), ignore.case = TRUE)
Old <- grep("old", unlist(Murine_HSC_Data[1,]), ignore.case = TRUE)

STHSC <- grep("STHSC", unlist(Murine_HSC_Data[1,]), ignore.case = TRUE)
LTHSC <- grep("LTHSC", unlist(Murine_HSC_Data[1,]), ignore.case = TRUE)
MPP <- grep("MPP", unlist(Murine_HSC_Data[1,]), ignore.case = TRUE)



AllData <- data.matrix(Murine_HSC_Data[-1, c(Old, Young)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

Categories <- rep(NA, max(c(Old, Young)))
Categories[intersect(MPP, Old)] <- "Old MPP"
Categories[intersect(MPP, Young)] <- "Young MPP"
Categories[intersect(LTHSC, Old)] <- "Old LTHSC"
Categories[intersect(LTHSC, Young)] <- "Young LTHSC"
Categories[intersect(STHSC, Old)] <- "Old STHSC"
Categories[intersect(STHSC, Young)] <- "Young STHSC"
Categories <- Categories[!is.na(Categories)]

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories


KowalczykInfo.DBA.All <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Lasso', DistillThr = .5,
                                               IgnoreTail = FALSE, Log = TRUE, StartQuant = .4, Title = "Kowalczyk et al (All)", PlanVarLimit = .85, PlanVarLimitIC = .95)

KowalczykInfo.DBA.All$WitPT
KowalczykInfo.DBA.All$FreePT








AllData <- data.matrix(Murine_HSC_Data[-1, LTHSC])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

Categories <- rep(NA, max(LTHSC))
# Categories[intersect(MPP, Old)] <- "Old MPP"
# Categories[intersect(MPP, Young)] <- "Young MPP"
Categories[intersect(LTHSC, Old)] <- "Old LTHSC"
Categories[intersect(LTHSC, Young)] <- "Young LTHSC"
# Categories[intersect(STHSC, Old)] <- "Old STHSC"
# Categories[intersect(STHSC, Young)] <- "Young STHSC"
Categories <- Categories[!is.na(Categories)]

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories


KowalczykInfo.DBA.LTHSC <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Lasso', DistillThr = .5,
                                                 IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al (All)", PlanVarLimit = .85, PlanVarLimitIC = .85)

KowalczykInfo.DBA.LTHSC$WitPT
KowalczykInfo.DBA.LTHSC$FreePT











AllData <- data.matrix(Murine_HSC_Data[-1, STHSC])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

Categories <- rep(NA, max(STHSC))
# Categories[intersect(MPP, Old)] <- "Old MPP"
# Categories[intersect(MPP, Young)] <- "Young MPP"
# Categories[intersect(LTHSC, Old)] <- "Old LTHSC"
# Categories[intersect(LTHSC, Young)] <- "Young LTHSC"
Categories[intersect(STHSC, Old)] <- "Old STHSC"
Categories[intersect(STHSC, Young)] <- "Young STHSC"
Categories <- Categories[!is.na(Categories)]

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories


KowalczykInfo.DBA.STHSC <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Lasso', DistillThr = .5,
                                                 IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al (All)", PlanVarLimit = .85, PlanVarLimitIC = .85)

KowalczykInfo.DBA.STHSC$WitPT
KowalczykInfo.DBA.STHSC$FreePT








AllData <- data.matrix(Murine_HSC_Data[-1, MPP])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

Categories <- rep(NA, max(MPP))
Categories[intersect(MPP, Old)] <- "Old MPP"
Categories[intersect(MPP, Young)] <- "Young MPP"
# Categories[intersect(LTHSC, Old)] <- "Old LTHSC"
# Categories[intersect(LTHSC, Young)] <- "Young LTHSC"
# Categories[intersect(STHSC, Old)] <- "Old STHSC"
# Categories[intersect(STHSC, Young)] <- "Young STHSC"
Categories <- Categories[!is.na(Categories)]

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories


KowalczykInfo.DBA.MPP <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Lasso', DistillThr = .5,
                                               IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Kowalczyk et al (All)", PlanVarLimit = .85, PlanVarLimitIC = .85)

KowalczykInfo.DBA.MPP$WitPT
KowalczykInfo.DBA.MPP$FreePT












































# Klein et al ----------------------------------------------------------------

# BaseDir <- "~/Google Drive/Datasets/Klein et al - Murine Embryonic Stem Cells/"
# 
# Mouse_ES_Data_d0 <- read_csv(paste(BaseDir, "GSM1599494_ES_d0_main.csv.bz2", sep=''), col_names = FALSE)
# Mouse_ES_Data_d0_R1 <- read_csv(paste(BaseDir, "GSM1599495_ES_d0_biorep_techrep1.csv.bz2", sep=''), col_names = FALSE)
# Mouse_ES_Data_d0_R2 <- read_csv(paste(BaseDir, "GSM1599496_ES_d0_biorep_techrep2.csv.bz2", sep =''), col_names = FALSE)
# 
# Mouse_ES_Data_d2 <- read_csv(paste(BaseDir, "GSM1599497_ES_d2_LIFminus.csv.bz2", sep=''), col_names = FALSE)
# Mouse_ES_Data_d4 <- read_csv(paste(BaseDir, "GSM1599498_ES_d4_LIFminus.csv.bz2", sep=''), col_names = FALSE)
# Mouse_ES_Data_d7 <- read_csv(paste(BaseDir, "GSM1599499_ES_d7_LIFminus.csv.bz2", sep=''), col_names = FALSE)
# 
# Genes <- as.character(Mouse_ES_Data_d0$X1)
# Mouse_ES_Data_d0 <- data.frame(Mouse_ES_Data_d0[, -1])
# rownames(Mouse_ES_Data_d0) <- Genes
# 
# Genes <- as.character(Mouse_ES_Data_d0_R1$X1)
# Mouse_ES_Data_d0_R1 <- data.frame(Mouse_ES_Data_d0_R1[, -1])
# rownames(Mouse_ES_Data_d0_R1) <- Genes
# 
# Genes <- as.character(Mouse_ES_Data_d0_R2$X1)
# Mouse_ES_Data_d0_R2 <- data.frame(Mouse_ES_Data_d0_R2[, -1])
# rownames(Mouse_ES_Data_d0_R2) <- Genes
# 
# Genes <- as.character(Mouse_ES_Data_d2$X1)
# Mouse_ES_Data_d2 <- data.frame(Mouse_ES_Data_d2[, -1])
# rownames(Mouse_ES_Data_d2) <- Genes
# 
# Genes <- as.character(Mouse_ES_Data_d4$X1)
# Mouse_ES_Data_d4 <- data.frame(Mouse_ES_Data_d4[, -1])
# rownames(Mouse_ES_Data_d4) <- Genes
# 
# Genes <- as.character(Mouse_ES_Data_d7$X1)
# Mouse_ES_Data_d7 <- data.frame(Mouse_ES_Data_d7[, -1])
# rownames(Mouse_ES_Data_d7) <- Genes
# 
# 
# 
# GeneExprMat <- data.matrix(Mouse_ES_Data_d0)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- rep("d0", ncol(Mouse_ES_Data_d0))
# 
# FilteredGenes.Klein.d0.Circle <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = "VarPC", DistillThr = .5, FastReduce = FALSE,
#                                               ExtMode = 2, StopCrit = .98, OutThr = 3, nNodes = 40, VarThr = .99, Categories = Categories,
#                                               GraphType = 'Circle', PlanVarLimit = .90, PlanVarLimitIC = .95, ForceLasso = FALSE,
#                                               LassoCircInit = 10, MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 25, DipPVThr = 1e-4,
#                                               PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
# 
# 
# 
# GeneExprMat <- data.matrix(Mouse_ES_Data_d0_R1)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- rep("d0_r1", ncol(Mouse_ES_Data_d0_R1))
# 
# FilteredGenes.Klein.d0_r1.Circle <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = "VarPC", DistillThr = .5, FastReduce = FALSE,
#                                                  ExtMode = 2, StopCrit = .98, OutThr = 3, nNodes = 40, VarThr = .99, Categories = Categories,
#                                                  GraphType = 'Circle', PlanVarLimit = .90, PlanVarLimitIC = .95, ForceLasso = FALSE,
#                                                  LassoCircInit = 10, MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 25, DipPVThr = 1e-4,
#                                                  PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
# 
# 
# 
# 
# GeneExprMat <- data.matrix(Mouse_ES_Data_d0_R2)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- rep("d0_r2", ncol(Mouse_ES_Data_d0_R2))
# 
# FilteredGenes.Klein.d0_r2.Circle <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = "VarPC", DistillThr = .5, FastReduce = FALSE,
#                                                     ExtMode = 2, StopCrit = .98, OutThr = 3, nNodes = 40, VarThr = .99, Categories = Categories,
#                                                     GraphType = 'Circle', PlanVarLimit = .90, PlanVarLimitIC = .95, ForceLasso = FALSE,
#                                                     LassoCircInit = 10, MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 25, DipPVThr = 1e-4,
#                                                     PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
# 
# 
# 
# 
# 
# 
# GeneExprMat <- data.matrix(Mouse_ES_Data_d2)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- rep("d0_r2", ncol(Mouse_ES_Data_d2))
# 
# FilteredGenes.Klein.d2.Circle <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = "VarPC", DistillThr = .5, FastReduce = FALSE,
#                                                     ExtMode = 2, StopCrit = .98, OutThr = 3, nNodes = 40, VarThr = .99, Categories = Categories,
#                                                     GraphType = 'Circle', PlanVarLimit = .90, PlanVarLimitIC = .95, ForceLasso = FALSE,
#                                                     LassoCircInit = 10, MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 25, DipPVThr = 1e-4,
#                                                     PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
# 
# 
# 
# 
# GeneExprMat <- data.matrix(Mouse_ES_Data_d4)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- rep("d0_r2", ncol(Mouse_ES_Data_d4))
# 
# FilteredGenes.Klein.d4.Circle <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = "VarPC", DistillThr = .5, FastReduce = FALSE,
#                                                  ExtMode = 2, StopCrit = .98, OutThr = 3, nNodes = 40, VarThr = .99, Categories = Categories,
#                                                  GraphType = 'Circle', PlanVarLimit = .90, PlanVarLimitIC = .95, ForceLasso = FALSE,
#                                                  LassoCircInit = 10, MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 25, DipPVThr = 1e-4,
#                                                  PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
# 
# 
# 
# GeneExprMat <- data.matrix(Mouse_ES_Data_d7)
# StartSet <- MouseGenes_GOCellCycle
# Categories <- rep("d0_r2", ncol(Mouse_ES_Data_d7))
# 
# FilteredGenes.Klein.d7.Circle <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = "VarPC", DistillThr = .5, FastReduce = FALSE,
#                                                  ExtMode = 2, StopCrit = .98, OutThr = 3, nNodes = 40, VarThr = .99, Categories = Categories,
#                                                  GraphType = 'Circle', PlanVarLimit = .90, PlanVarLimitIC = .95, ForceLasso = FALSE,
#                                                  LassoCircInit = 10, MinBranDiff = 2, Log = TRUE, Filter = TRUE, MinProlCells = 25, DipPVThr = 1e-4,
#                                                  PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
# 
# 
# rm(Mouse_ES_Data_d0, Mouse_ES_Data_d0_R1, Mouse_ES_Data_d0_R2, Mouse_ES_Data_d2, Mouse_ES_Data_d4, Mouse_ES_Data_d7)






# Llorens-Bobadilla et al ----------------------------------------------------------------

BaseDir <- "~/Google Drive/Datasets/Llorens-Bobadilla et al. - Murine Adult Neural Stem Cells/"

GSE67833_Gene_expression_matrix_csv <- read_csv("~/Google Drive/Datasets/Llorens-Bobadilla et al. - Murine Adult Neural Stem Cells/GSE67833_Gene_expression_matrix.csv.gz")

Murine_ANSC_Data1 <- read_csv(paste(BaseDir, "GSE67833_Gene_expression_matrix.csv.gz", sep= ''))
Murine_ANSC_Data2 <- read_csv(paste(BaseDir, "GSE67833_Gene_expression_matrix_GSM1684656-704.csv.gz", sep= ''))

GSE67833_series_matrix_txt <- read_csv(paste(BaseDir, "GSE67833_series_matrix_proc.csv", sep = ''), trim_ws = TRUE)

SampleInfo_1 <- data.frame(GSE67833_series_matrix_txt[,1:224])
SampleInfo_2 <- data.frame(GSE67833_series_matrix_txt[,c(1,225:273)])

table(unlist(SampleInfo_1[7,-1]))
table(unlist(SampleInfo_1[9,-1]))

table(unlist(SampleInfo_2[7,-1]))
table(unlist(SampleInfo_2[14,-1]))


Names_1 <- colnames(SampleInfo_1)[-1]
Origin_2 <- unlist(SampleInfo_1[7,-1])
Condition_1 <- unlist(SampleInfo_1[9,-1])

Type_1 = rep("Other", length(Names_1))
Type_1[grep("PSA", Names_1)] = "PSA"
Type_1[grep("GP", Names_1)] = "GP"

table(Type_1, Condition_1)

Names_2 <- unlist(lapply(strsplit(colnames(SampleInfo_2)[-1], split = ".", fixed = TRUE), "[[", 3))
Origin_2 <- unlist(lapply(strsplit(unlist(SampleInfo_2[7,-1]), "_"), paste, collapse=' '))
Condition_2 <- unlist(SampleInfo_2[14,-1])

colnames(Murine_ANSC_Data1)[1] <- "Gene"
colnames(Murine_ANSC_Data2)[1] <- "Gene"

GName1 <- Murine_ANSC_Data1$Gene
GName2 <- Murine_ANSC_Data2$Gene

Murine_ANSC_Data1 <- data.matrix(Murine_ANSC_Data1[,-1])
Murine_ANSC_Data2 <- data.matrix(Murine_ANSC_Data2[,-1])

rownames(Murine_ANSC_Data1) <- GName1
rownames(Murine_ANSC_Data2) <- GName2


ncol(Murine_ANSC_Data1)
nrow(Murine_ANSC_Data1)
ncol(Murine_ANSC_Data2)
nrow(Murine_ANSC_Data2)


GeneExprMat <- data.matrix(Murine_ANSC_Data1)
StartSet <- ConvertNames(SourceOrganism = "mouse", TargetOrganism = "mouse", Genes = MouseGenes_GOCellCycle, SourceTypes = "Names", TargetTypes = "Ensembl")
Categories <- factor(Condition_1)

FilteredGenes.Llor.Circle <- PlotsAndDistill.Mouse(GeneExprMat = GeneExprMat, StartSet = StartSet, Categories = Categories, Topology = 'Lasso', DistillThr = .5,
                                                   IgnoreTail = FALSE, Log = TRUE, StartQuant = .5, Title = "Llorens-Bobadilla et al (All)", PlanVarLimit = .85, PlanVarLimitIC = .9)

FilteredGenes.Llor.Circle$WitPT
FilteredGenes.Llor.Circle$FreePT
















# Shin et al ----------------------------------------------------------------


# Treutlein et al ----------------------------------------------------------------


# Usoskin et al ----------------------------------------------------------------




















# HUMAN









