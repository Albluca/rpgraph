# 
# 
# # Trapnell
# 
# 
# 
# library(readr)
# 
# BaseDir <- "~/Google Drive/Datasets/Trapnell et al - Human skeletal muscle myoblasts/"
# 
# Data_Header <- read_delim(paste(BaseDir, "GSE52529_fpkm_matrix.txt.gz", sep = ''), "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, n_max = 1)
# Data_Myo <- read_delim(paste(BaseDir, "GSE52529_fpkm_matrix.txt.gz", sep = ''), "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, skip = 1)
# 
# Data_Header_2 <- read_delim(paste(BaseDir, "GSE52529_truseq_fpkm_matrix.txt.gz", sep = ''), "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, n_max = 1)
# Data_Myo_2 <- read_delim(paste(BaseDir, "GSE52529_truseq_fpkm_matrix.txt.gz", sep = ''), "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE, skip = 1)
# 
# Sample_InfoHeader <-  read_delim(paste(BaseDir, "GSE52529-GPL11154_series_matrix.txt.gz", sep = ''), "\t", escape_double = FALSE, n_max = 38,
#                                  col_names = FALSE, trim_ws = TRUE)
# Sample_Info <-  read_delim(paste(BaseDir, "GSE52529-GPL11154_series_matrix.txt.gz", sep = ''), "\t", escape_double = FALSE, skip = 38, 
#                            col_names = FALSE, trim_ws = TRUE, n_max = 47)
# 
# 
# Sample_InfoHeader_2 <-  read_delim(paste(BaseDir, "GSE52529-GPL16791_series_matrix.txt.gz", sep = ''), "\t", escape_double = FALSE, n_max = 38,
#                                    col_names = FALSE, trim_ws = TRUE)
# Sample_Info_2 <-  read_delim(paste(BaseDir, "GSE52529-GPL16791_series_matrix.txt.gz", sep = ''), "\t", escape_double = FALSE, skip = 38, 
#                              col_names = FALSE, trim_ws = TRUE, n_max = 47)
# 
# 
# 
# table(unlist(Sample_Info[11,-1]))
# table(unlist(Sample_Info_2[11,-1]))
# 
# all(Data_Myo[,1] == Data_Myo_2[,1])
# 
# # Removing the decimal part. Not sure why this is necessary 
# 
# # Genes <- sapply(strsplit(unlist(Data_Myo[,1]), ".", TRUE), "[[", 1)
# 
# # library(biomaRt)
# 
# AllData <- as.data.frame(cbind(Data_Myo[,-1], Data_Myo_2[,-1]))
# 
# Names <- c(unlist(Data_Header[1,]), unlist(Data_Header_2[1,]))
# 
# colnames(AllData) <- Names
# rownames(AllData) <- unlist(Data_Myo[,1])
# 
# ToRead <- grep("_CT_", Names)
# 
# Names <- Names[ToRead]
# AllData <- AllData[,ToRead]
# 
# TimeInfo <- sapply(strsplit(Names, "_"), "[[", 1)
# table(TimeInfo)
# 
# 
# 
# library(rpgraph)
# 
# AllData <- AllData[(rowSums(AllData==0)/ncol(AllData))<.75, ]
# 
# OutExpr <- scater::isOutlier(rowSums(AllData>0), nmads = 5)
# OutCount <- scater::isOutlier(rowSums(AllData), nmads = 5)
# 
# AllData <- AllData[!(OutCount | OutExpr), ]
# 
# dim(AllData)
# 
# PCAData <- prcomp(t(AllData), retx = TRUE, scale. = FALSE, center = FALSE)
# 
# plot(PCAData$x[,1:2])
# 
# Centroid <- colMeans(PCAData$x)
# 
# Dists <- as.matrix(dist(rbind(Centroid, PCAData$x)))
# DistFromCent <- Dists[1,-1]
# 
# PCAFil <- scater::isOutlier(DistFromCent, nmads = 5)
# 
# 
# plot(PCAData$x[!PCAFil,1:2])
# 
# 
# 
# BasicTreeData <- computeElasticPrincipalGraph(Data = PCAData$x[!PCAFil,], NumNodes = 60, Method = 'DefaultPrincipalTreeConfiguration', NodeStep = 1)
# 
# plotPieNet(Results = BasicTreeData[[length(BasicTreeData)]], Data = PCAData$x[!PCAFil,], Categories = factor(TimeInfo[!PCAFil]))
# 
# 


# Schlitzer




PlotsAndDistill.DC <- function(GeneExprMat, StartSet, Categories, Topology = "Line", PlanVarLimit = .9, PlanVarLimitIC = 92,
                                  DistillThr = .7, Log = TRUE, StartQuant = .5, Title = '', IgnoreTail = FALSE) {
  
  PrGraph.Initial <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = StartSet, OutThr = 3, nNodes = 20,
                                       VarThr = .99, Categories = Categories, GraphType = Topology, PlanVarLimit = PlanVarLimit,
                                       PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, InitStructNodes = 10,
                                       MinBranDiff = 2, Log = Log, Filter = TRUE, MinProlCells = 200, DipPVThr = 1e-4,
                                       PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
  
  
  
  # DataSet = GeneExprMat
  # GeneSet = StartSet
  # OutThr = 3
  # nNodes = 40
  # VarThr = .99
  # Categories = Categories
  # GraphType = Topology
  # PlanVarLimit = PlanVarLimit
  # PlanVarLimitIC = PlanVarLimitIC
  # ForceLasso = FALSE
  # LassoCircInit = 20
  # MinBranDiff = 2
  # Log = Log
  # Filter = TRUE
  # MinProlCells = 20
  # DipPVThr = 1e-4
  # PCACenter = FALSE
  # PCAProjCenter = TRUE
  # PlotIntermediate = FALSE
  
  
  DistilledGenes <- ConvergeOnGenes(ExpData = GeneExprMat, StartGeneSet = StartSet, Mode = "VarPC", DistillThr = DistillThr, FastReduce = FALSE,
                                    ExtMode = 2, StopCrit = .95, OutThr = 3, nNodes = 40, VarThr = .99, Categories = Categories,
                                    GraphType = Topology, PlanVarLimit = PlanVarLimit, PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, 
                                    MinBranDiff = 2, Log = Log, Filter = TRUE, MinProlCells = 200, DipPVThr = 1e-4,
                                    PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE,
                                    IgnoreTail = IgnoreTail, StartQuant = StartQuant)
  
  PrGraph.Final <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = DistilledGenes, OutThr = 3, nNodes = 40,
                                     VarThr = .99, Categories = Categories, GraphType = Topology, PlanVarLimit = PlanVarLimit,
                                     PlanVarLimitIC = PlanVarLimitIC, ForceLasso = FALSE, InitStructNodes = 20,
                                     MinBranDiff = 2, Log = Log, Filter = TRUE, MinProlCells = 200, DipPVThr = 1e-4,
                                     PCACenter = FALSE, PCAProjCenter = TRUE, PlotIntermediate = FALSE)
  
  
  
  
  ProjectOnPrincipalGraph(Nodes = PrGraph.Initial$PrinGraph$Nodes, Edges = PrGraph.Initial$PrinGraph$Edges, Points = PrGraph.Initial$Data,
                          UsedPoints = NULL, Categories = PrGraph.Initial$Categories, Title=paste(Title, "(All Genes)"),
                          PCACenter = TRUE, ShowFitted = FALSE)
  
  ProjectOnPrincipalGraph(Nodes = PrGraph.Final$PrinGraph$Nodes, Edges = PrGraph.Final$PrinGraph$Edges, Points = PrGraph.Final$Data,
                          UsedPoints = NULL, Categories = PrGraph.Final$Categories, Title=paste(Title, "(Filtered)"),
                          PCACenter = TRUE, ShowFitted = FALSE)
  
  return(list(Genes = DistilledGenes, FinalStruct = PrGraph.Final))
  
}






library(GEOquery)
gse <- getGEO("GSE60783", GSEMatrix = TRUE)

filePaths = getGEOSuppFiles("GSE60783")

library(readr)
Content <- untar(row.names(filePaths)[1], list = TRUE)
untar(row.names(filePaths)[1], list = FALSE)

MatData <- NULL
UsedFiles <- NULL

for(i in 1:length(Content)){
  if(length(grep(pattern = "readcount", x = Content[i], ignore.case = TRUE)) == 0){
    next()
  } else {
    UsedFiles <- c(UsedFiles, i)
  }
  
  Exp <- read_delim(Content[i], "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
  
  if(is.null(MatData)){
    MatData <- cbind(unlist(Exp[,1]), unlist(Exp[,2]))
  } else {
    if(any(MatData[,1] != unlist(Exp[,1]))){
      stop("Incompatible samples")
    }
    MatData <- cbind(MatData, unlist(Exp[,2]))
  }
  
  file.remove(Content[i])
  
}

SplitPath <- unlist(strsplit(rownames(filePaths)[1], "/"))
SplitPath <- SplitPath[-length(SplitPath)]
unlink(x = paste(SplitPath, collapse = "/"), recursive = TRUE)

Genes <- MatData[,1]
MatData <- data.matrix(data.frame(MatData[,-1]))
rownames(MatData) <- Genes
colnames(MatData) <- unlist(lapply(strsplit(Content[UsedFiles], "_"), "[[", 1))


SampIDs <- as.character(gse$`GSE60783-GPL13112_series_matrix.txt.gz`[[2]])

CellPop <- as.character(gse$`GSE60783-GPL13112_series_matrix.txt.gz`[[12]])

CellPop <- unlist(lapply(strsplit(CellPop, ": "), "[[", 2))
CellPop <- unlist(lapply(strsplit(CellPop, "_sorted"), "[[", 1))
CellPop <- factor(CellPop, levels = c("MDP", "CDP", "PreDC"))


CellPop[match(colnames(MatData), SampIDs)]


CellPop










library(rpgraph)




ProcessedData <- PlotsAndDistill.DC(GeneExprMat = MatData, StartSet = rownames(MatData), Categories = CellPop,
                                    Topology = "Line", PlanVarLimit = .85, PlanVarLimitIC = 9,
                               DistillThr = .4, Log = TRUE, StartQuant = .5, Title = '')
  


PlotOnStages(Structure = "Line", TaxonList = ProcessedData$FinalStruct$TaxonList[[length(ProcessedData$FinalStruct$TaxonList)]],
             Categories = ProcessedData$FinalStruct$Categories, nGenes = 50,
             PrinGraph = ProcessedData$FinalStruct$PrinGraph,
             Net = ProcessedData$FinalStruct$Net[[length(ProcessedData$FinalStruct$Net)]],
             SelThr = .35, ComputeOverlaps = TRUE, ExpData = ProcessedData$FinalStruct$FiltExp,
             RotatioMatrix = ProcessedData$FinalStruct$PCAData$rotation[,1:ProcessedData$FinalStruct$nDims],
             PointProjections = ProcessedData$FinalStruct$ProjPoints[[length(ProcessedData$FinalStruct$ProjPoints)]])

ColCat <- c("green", "red", "blue")
names(ColCat) <- c("MDP", "CDP", "PreDC")


plotPieNet(Results =  ProcessedData$FinalStruct$IntGrahs[[length(ProcessedData$FinalStruct$IntGrahs)]],
                            Data = ProcessedData$FinalStruct$Data, NodeSizeMult = 4,
                            Categories = ProcessedData$FinalStruct$Categories, PlotNet = TRUE,
                            Graph = ProcessedData$FinalStruct$Net[[length(ProcessedData$FinalStruct$Net)]],
                            TaxonList = ProcessedData$FinalStruct$TaxonList[[length(ProcessedData$FinalStruct$TaxonList)]],
                            LayOut = 'circle_line', Main = "Pincipal graph", ColCat = ColCat)

legend(x = "center", legend = names(ColCat), fill = ColCat)

# 
# 
# 
# Structure = "Line"
# TaxonList = ProcessedData$FinalStruct$TaxonList[[length(ProcessedData$FinalStruct$TaxonList)]]
# Categories = ProcessedData$FinalStruct$Categories
# nGenes = 20
# PrinGraph = ProcessedData$FinalStruct$PrinGraph
# Net = ProcessedData$FinalStruct$Net[[length(ProcessedData$FinalStruct$Net)]]
# SelThr = .35
# ComputeOverlaps = TRUE
# ExpData = ProcessedData$FinalStruct$FiltExp
# RotatioMatrix = ProcessedData$FinalStruct$PCAData$rotation[,1:ProcessedData$FinalStruct$nDims]
# PointProjections = ProcessedData$FinalStruct$ProjPoints[[length(ProcessedData$FinalStruct$ProjPoints)]]
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
# AllData <- log10(MatData + 1)
# 
# # plot(apply(AllData, 1, sd), apply(AllData, 1, mean))
# 
# VarVect <- apply(AllData, 1, var)
# 
# AllData <- AllData[VarVect > quantile(VarVect, .25), ]
# 
# # AllData <- AllData[(rowSums(AllData==0)/ncol(AllData))<.75, ]
# 
# # OutExpr <- scater::isOutlier(rowSums(AllData>0), nmads = 5)
# # OutCount <- scater::isOutlier(rowSums(AllData), nmads = 5)
# 
# # AllData <- AllData[!(OutCount | OutExpr), ]
# 
# 
# 
# dim(AllData)
# 
# PCAData <- prcomp(t(AllData), retx = TRUE, scale. = FALSE, center = FALSE)
# 
# plot(PCAData$x[,1:2])
# 
# Centroid <- colMeans(PCAData$x)
# 
# Dists <- as.matrix(dist(rbind(Centroid, PCAData$x)))
# DistFromCent <- Dists[1,-1]
# 
# PCAFil <- scater::isOutlier(DistFromCent, nmads = 5)
# 
# 
# plot(PCAData$x[!PCAFil,1:2])
# 
# 
# 
# 
# 
# 
# 
# BasicTreeData <- computeElasticPrincipalGraph(Data = PCAData$x[!PCAFil,], NumNodes = 15, Method = 'CurveConfiguration', NodeStep = 1)
# 
# Net <- ConstructGraph(Results = BasicTreeData[[length(BasicTreeData)]], DirectionMat = NULL, Thr = 0.05)
# 
# TaxonList <- getTaxonMap(Results = BasicTreeData[[length(BasicTreeData)]], Data = PCAData$x[!PCAFil,])
# 
# PieINFO <- plotPieNet(Results = BasicTreeData[[length(BasicTreeData)]], Data = PCAData$x[!PCAFil,], Categories = CellPop[!PCAFil],
#                       NodeSizeMult = 2, Graph = Net, TaxonList = TaxonList)
# PieINFO$ColInfo
# 
# legend("bottomleft", legend = names(PieINFO$ColInfo), fill = PieINFO$ColInfo)
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
# ProjPoints <- projectPoints(Results = BasicTreeData[[length(BasicTreeData)]], Data = PCAData$x[!PCAFil,],
#                             TaxonList = TaxonList, UseR = TRUE)
# 
# Path <- igraph::graph.get.subisomorphisms.vf2(graph2 = igraph::graph.ring(n = 15, directed = FALSE, circular = FALSE), graph1 = PieINFO$Net)
# 
# 
# TaxVect <- rep(NA, max(unlist(TaxonList), na.rm = TRUE))
# for(i in 1:length(TaxonList)){
#   TaxVect[TaxonList[[i]]] <- i 
# }
# 
# NumPath <- unlist(lapply(strsplit(Path[[1]]$name, "V_"), "[[", 2))
# 
# Reordered <- TaxVect
# for(j in 1:length(NumPath)){
#   Reordered[TaxVect == NumPath[j]] <- j
# }
# 
# 
# 
# boxplot(Reordered ~ factor(CellPop[!PCAFil]))
# 
# 
# Empty <- which(unlist(lapply(lapply(TaxonList, is.na), any)))
# 
# TB <- table(CellPop[!PCAFil], factor(TaxVect, levels = paste(1:15)))
# colnames(TB)
# 
# 
# NodeOnGenes <- t(BasicTreeData[[length(BasicTreeData)]]$Nodes %*% t(PCAData$rotation))
# 
# 
# SelPath <- NumPath
# # SelPath <- SelPath[-length(SelPath)]
# 
# Idx <- 7
# 
# library(ggplot2)
# 
# ggplot(data = data.frame(x=1:ncol(NodeOnGenes), y=NodeOnGenes[Idx,as.numeric(SelPath)]), mapping = aes(x = x, y = y)) +
#   geom_point() + geom_line() + labs(x = "Pseudotime", y="Gene expression (log10 Pseudocount)", title = rownames(NodeOnGenes)[Idx]) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# 
# 
# 
# 
# TB <- TB[, NumPath]
# 
# barplot(t(t(TB)/colSums(TB)), col = c("red", "green", "blue"), beside = TRUE, las = 2)
# 
# 
# 
# 
# PercMat <- t(t(TB)/colSums(TB))
# PercMat[is.na(PercMat)] <- 0
# BinPercMat <- (PercMat > .4)
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
# 
# for(i in 2:(ncol(BinPercMat)-1)){
#   for(j in 1:nrow(BinPercMat)){
#     if(BinPercMat[j,i-1] & BinPercMat[j,i+1]){
#       BinPercMat[j,i] <- TRUE
#     }
#   }
# }
# 
# 
# 
# for(i in 2:(ncol(BinPercMat)-1)){
#   for(j in 1:nrow(BinPercMat)){
#     if(!BinPercMat[j,i-1] & !BinPercMat[j,i+1]){
#       BinPercMat[j,i] <- FALSE
#     }
#   }
# }
# 
# 
# 
# rownames(BinPercMat) <- c("MDP", "CDP", "PreDC") 
# 
# 
# 
# P12 <- BinPercMat[1,] & BinPercMat[2,]
# 
# P23 <- BinPercMat[2,] & BinPercMat[3,]
# 
# BinPercMatExt <- rbind(BinPercMat[1,] & !P12, P12,
#                        BinPercMat[2,] & !P12 & !P23, P23,
#                        BinPercMat[3,])*1
# rownames(BinPercMatExt) <- c("MDP", "MDP/CDP", "CDP", "CDP/PreDC", "PreDC")
# 
# 
# 
# 
# NodeOnGenes <- t(BasicTreeData[[length(BasicTreeData)]]$Nodes %*% t(PCAData$rotation))
# 
# OrderedPoints <- OrderOnPath(PrinGraph = BasicTreeData[[length(BasicTreeData)]], Path = NumPath, PointProjections = ProjPoints)
# 
# 
# Idxs <- apply(BinPercMatExt==1, 1, which)
# 
# Bond <- lapply(Idxs[lapply(Idxs, length) > 1], range)
# 
# 
# for(Idx in sample(1:nrow(NodeOnGenes), 10)){
#   
#   p <- ggplot(data = data.frame(x=cumsum(OrderedPoints$PathLen),
#                                 y=NodeOnGenes[Idx,as.numeric(NumPath)]),
#               mapping = aes(x = x, y = y, color="PC")) + labs(x = "Pseudotime", y="Gene expression (log10 Pseudocount)", title = rownames(NodeOnGenes)[Idx]) +
#     theme(plot.title = element_text(hjust = 0.5))
#   
#   RangeDF <- data.frame(t(as.data.frame(lapply(Idxs, range))))
#   colnames(RangeDF) <- c("Min", "Max")
#   RangeDF[!is.finite(RangeDF$Min),] <- NA
#   RangeDF <- cbind(rownames(RangeDF), RangeDF)
#   colnames(RangeDF) <- c("Stage", "Min", "Max")
#   
#   RangeDF$Stage <- factor(as.character(RangeDF$Stage), levels = c("MDP", "MDP.CDP", "CDP", "CDP.PreDC", "PreDC"))
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
#     geom_point(data = data.frame(x=OrderedPoints$PositionOnPath, y=AllData[Idx, !PCAFil]),
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
# 
# 
# 
# 
# 
# 
# # 
# # 
# # #  Setty et al
# # 
# # library(readr)
# # Setty_ColNames <- read_csv("~/Google Drive/Datasets/Setty et al/mouse_marrow_scrnaseq.csv.zip", n_max = 1, col_names = FALSE)
# # Setty_Data <- read_csv("~/Google Drive/Datasets/Setty et al/mouse_marrow_scrnaseq.csv.zip", skip = 1, col_names = FALSE)
# # 
# # CellNames <- unlist(Setty_Data[,1])
# # CellInfo <- data.matrix(Setty_Data[,(ncol(Setty_Data)-3):ncol(Setty_Data)])
# # colnames(CellInfo) <- unlist(Setty_ColNames[c((ncol(Setty_Data)-3):ncol(Setty_Data))])
# # rownames(CellInfo) <- CellNames
# # 
# # Setty_Data <- data.matrix(Setty_Data[,-c(1,(ncol(Setty_Data)-3):ncol(Setty_Data))])
# # colnames(Setty_Data) <- unlist(Setty_ColNames[-c(1,(ncol(Setty_Data)-3):ncol(Setty_Data))])
# # rownames(Setty_Data) <- CellNames
# # 
# # 
# # 
# # 
# # plot(apply(Setty_Data, 2, var), apply(Setty_Data, 2, mean))
# # 
# # 
# # 
# # dim(Setty_Data)
# # 
# # PCAData <- prcomp(Setty_Data, retx = TRUE, scale. = FALSE, center = FALSE)
# # 
# # plot(PCAData$x[,1:2])
# # 
# # Centroid <- colMeans(PCAData$x)
# # 
# # Dists <- as.matrix(dist(rbind(Centroid, PCAData$x)))
# # DistFromCent <- Dists[1,-1]
# # 
# # PCAFil <- scater::isOutlier(DistFromCent, nmads = 5)
# # sum(PCAFil)
# # 
# # plot(PCAData$x[!PCAFil,1:2])
# # 
# # library(rpgraph)
# # 
# # 
# # BasicTreeData <- computeElasticPrincipalGraph(Data = PCAData$x[!PCAFil,], NumNodes = 12, Method = 'DefaultPrincipalTreeConfiguration', NodeStep = 1)
# # 
# # PieINFO <- plotPieNet(Results = BasicTreeData[[length(BasicTreeData)]], Data = PCAData$x[!PCAFil,], Categories = factor(CellPop[!PCAFil]))
# # PieINFO$ColInfo
# # 
# # legend("bottomleft", legend = names(PieINFO$ColInfo), fill = PieINFO$ColInfo)
# # 
# # 
# # 
# # 
# # 
# # 
