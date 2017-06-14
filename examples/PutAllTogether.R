# load("/Users/newmac-luca/Desktop/RecentData_U.RData")
# 
# FinalStructure.Kowa <- ExpDataset.Kowa$ListProc[[length(ExpDataset.Kowa$ListProc)]]
# TaxonList.Kowa <- ExpDataset.Kowa$StartInfo$FinalStruct$TaxonList[[length(ExpDataset.Kowa$StartInfo$FinalStruct$TaxonList)]]
# ReOrd.Kowa <- match(colnames(FinalStructure.Kowa$CellExp), names(FinalStructure.Kowa$CellsPT))
# Genes.Kowa <- sort(unique(c(rownames(FinalStructure.Kowa$NodesExp), rownames(FinalStructure.Kowa$NodesExp))))
# 
# 
# FinalStructure.Sasa <- ExpDataset.Sasa$ListProc[[length(ExpDataset.Sasa$ListProc)]]
# TaxonList.Sasa <- ExpDataset.Sasa$StartInfo$FinalStruct$TaxonList[[length(ExpDataset.Sasa$StartInfo$FinalStruct$TaxonList)]]
# ReOrd.Sasa <- match(colnames(FinalStructure.Sasa$CellExp), names(FinalStructure.Sasa$CellsPT))
# Genes.Sasa <- sort(unique(c(rownames(FinalStructure.Sasa$NodesExp), rownames(FinalStructure.Sasa$NodesExp))))
# 
# 
# FinalStructure.Buett <- ExpDataset.Buett$ListProc[[length(ExpDataset.Buett$ListProc)]]
# TaxonList.Buett <- ExpDataset.Buett$StartInfo$FinalStruct$TaxonList[[length(ExpDataset.Buett$StartInfo$FinalStruct$TaxonList)]]
# ReOrd.Buett <- match(colnames(FinalStructure.Buett$CellExp), names(FinalStructure.Buett$CellsPT))
# Genes.Buett <- sort(unique(c(rownames(FinalStructure.Buett$NodesExp), rownames(FinalStructure.Buett$NodesExp))))
# 
# 
# tGMT <- rRoma::SelectFromMSIGdb("GO_CELL_CYCLE")
# MouseGenes_GOCellCycle <- tGMT[[which(unlist(lapply(tGMT, "[[", "Name"))=="GO_CELL_CYCLE")]]$Genes
# MouseGenes_GOCellCycle <- gsub("(\\b[a-z]{1})", "\\U\\1", tolower(MouseGenes_GOCellCycle), perl=TRUE)

# options("java.home"="/Library/Java/JavaVirtualMachines/jdk1.8.0_112.jdk/Contents/Home/jre")
# dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_112.jdk/Contents/Home/jre/lib/server/libjvm.dylib')

library(rpgraph)

Data.Sasa <- read_rds("~/Google Drive/Datasets/Sasagawa et al - Murine Stem Cells (EB5 cell line)/rPG_Last.rds")
Data.Kowa <- read_rds("~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.rds")
Data.Buet <- read_rds("~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/rPG_Last.rds")



TargetStruct <- Data.Kowa$Analysis$PGStructs[[length(Data.Kowa$Analysis$PGStructs)]]
Proc.Exp.Kowa <- PlotOnStages(Structure = "Circle",
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
                         OrderOnCat = TRUE,
                         SmoothPoints = 2, MinCellPerNode = 2)


TargetStruct <- Data.Buet$Analysis$PGStructs[[length(Data.Buet$Analysis$PGStructs)]]
Proc.Exp.Buet <- PlotOnStages(Structure = "Circle",
                              Categories = TargetStruct$Categories,
                              nGenes = 2,
                              TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                              PrinGraph = TargetStruct$PrinGraph,
                              Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                              SelThr = .35,
                              ComputeOverlaps = TRUE,
                              ExpData = TargetStruct$FiltExp,
                              RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                              PCACenter = TargetStruct$PCAData$center,
                              PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                              OrderOnCat = TRUE,
                              SmoothPoints = 1, MinCellPerNode = 2)


TargetStruct <- Data.Sasa$Analysis$PGStructs[[length(Data.Sasa$Analysis$PGStructs)]]
Proc.Exp.Sasa <- PlotOnStages(Structure = "Circle",
                              Categories = TargetStruct$Categories,
                              nGenes = 2,
                              TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                              PrinGraph = TargetStruct$PrinGraph,
                              Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                              SelThr = .35,
                              ComputeOverlaps = TRUE,
                              ExpData = TargetStruct$FiltExp,
                              RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                              PCACenter = TargetStruct$PCAData$center,
                              PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                              OrderOnCat = TRUE,
                              SmoothPoints = 1, MinCellPerNode = 1)





InputList <- list(list(Name = "Buettner et al", Expression = Data.Buet$ExpMat,
                       Categories = Data.Buet$Cats,
                       OrderedData = Proc.Exp.Buet,
                       PGStruct = Data.Buet$Analysis$PGStructs[[length(Data.Buet$Analysis$PGStructs)]],
                       TaxonList = Data.Buet$Analysis$PGStructs[[length(Data.Buet$Analysis$PGStructs)]]$TaxonList
                       ),
                  list(Name = "Kowalczyk et al", Expression = Data.Kowa$ExpMat,
                       Categories = Data.Kowa$Cats,
                       OrderedData = Proc.Exp.Kowa,
                       PGStruct = Data.Kowa$Analysis$PGStructs[[length(Data.Kowa$Analysis$PGStructs)]],
                       TaxonList = Data.Kowa$Analysis$PGStructs[[length(Data.Kowa$Analysis$PGStructs)]]$TaxonList
                  ),
                  list(Name = "Sasagawa et al", Expression = Data.Buet$ExpMat,
                       Categories = Data.Sasa$Cats,
                       OrderedData = Proc.Exp.Sasa,
                       PGStruct = Data.Sasa$Analysis$PGStructs[[length(Data.Sasa$Analysis$PGStructs)]],
                       TaxonList = Data.Sasa$Analysis$PGStructs[[length(Data.Sasa$Analysis$PGStructs)]]$TaxonList
                  )
                )



# InputList <- list(list(Name = "Buettner et al", Expression = FullExpData.Buett,
#                        Categories = FullCat.Buett, PrinCurveStruct = ExpDataset.Buett,
#                        TaxonList = TaxonList.Buett),
#                   list(Name = "Kowalczyk et al", Expression = FullExpData.Kowa,
#                        Categories = FullCat.Kowa, PrinCurveStruct = ExpDataset.Kowa,
#                        TaxonList = TaxonList.Kowa),
#                   list(Name = "Sasagawa et al", Expression = FullExpData.Sasa,
#                        Categories = FullCat.Sasa, PrinCurveStruct = ExpDataset.Sasa,
#                        TaxonList = TaxonList.Sasa))
# 
# rm(FullCat.Buett, FullCat.Kowa, FullCat.Sasa,
#    ExpDataset.Buett, ExpDataset.Kowa, ExpDataset.Sasa,
#    FullExpData.Buett, FullExpData.Kowa, FullExpData.Sasa)



CompareAcrossData(InputList)


# i <- 1
# 
# InputList[[i]]$TaxonList
# TaxVect <- rep(NA, length(InputList[[i]]$TaxonList))
# 
# for(j in 1:length(InputList[[i]]$TaxonList)){
#   TaxVect[InputList[[i]]$TaxonList[[j]]] <- j
# }
# 
# WorkStruct <- InputList[[i]]$PrinCurveStruct$ListProc[[length(InputList[[i]]$PrinCurveStruct$ListProc)]]
# SampleReord <- lapply(as.list(1:100), function(i){sample(TaxVect)})
# PvVect <- rep(NA, nrow(InputList[[i]]$Expression))
# 
# DoStuff <- function(gId) {
#   Exp <- unlist(InputList[[i]]$Expression[gId,names(WorkStruct$CellsPT)])
#   Base <- median(unlist(lapply(split(Exp, f = factor(TaxVect)), mad)))
# 
#   MedVect <- sapply(SampleReord, function(x){
#     median(unlist(lapply(split(Exp, f = factor(x)), mad)))
#   })
#   return(wilcox.test(MedVect - Base, alternative = "greater")$p.value)
# }
# 
# cl <- parallel::makeCluster(4, type = "FORK")
# parallel::clusterExport(cl, varlist = c("InputList", "WorkStruct", "TaxVect", "SampleReord", "i"), envir = environment())
# 
# AllPV.Buett <- parallel::parLapply(cl, as.list(1:nrow(InputList[[i]]$Expression)), DoStuff)
# 
# 
# AllPV.Buett <- unlist(AllPV.Buett)
# names(AllPV.Buett) <- rownames(InputList[[i]]$Expression)
# 
# 
# write_rds(AllPV.Buett, path = "AllPV.Buett.rds")
# 


# 
# 
# 
# 
# 
# 
# 
# 
# i <- 2
# 
# InputList[[i]]$TaxonList
# TaxVect <- rep(NA, length(InputList[[i]]$TaxonList))
# 
# for(j in 1:length(InputList[[i]]$TaxonList)){
#   TaxVect[InputList[[i]]$TaxonList[[j]]] <- j
# }
# 
# WorkStruct <- InputList[[i]]$PrinCurveStruct$ListProc[[length(InputList[[i]]$PrinCurveStruct$ListProc)]]
# SampleReord <- lapply(as.list(1:100), function(i){sample(TaxVect)})
# PvVect <- rep(NA, nrow(InputList[[i]]$Expression))
# 
# DoStuff <- function(gId) {
#   Exp <- unlist(InputList[[i]]$Expression[gId,names(WorkStruct$CellsPT)])
#   Base <- median(unlist(lapply(split(Exp, f = factor(TaxVect)), mad)))
#   
#   MedVect <- sapply(SampleReord, function(x){
#     median(unlist(lapply(split(Exp, f = factor(x)), mad)))
#   })
#   return(wilcox.test(MedVect - Base, alternative = "greater")$p.value)
# }
# 
# cl <- parallel::makeCluster(4, type = "FORK")
# parallel::clusterExport(cl, varlist = c("InputList", "WorkStruct", "TaxVect", "SampleReord", "i"), envir = environment())
# 
# AllPV.Kowa <- parallel::parLapply(cl, as.list(1:nrow(InputList[[i]]$Expression)), DoStuff)
# 
# 
# AllPV.Kowa <- unlist(AllPV.Kowa)
# names(AllPV.Kowa) <- rownames(InputList[[i]]$Expression)
# 
# 
# write_rds(AllPV.Kowa, path = "AllPV.Kowa.rds")
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
# i <- 3
# 
# InputList[[i]]$TaxonList
# TaxVect <- rep(NA, length(InputList[[i]]$TaxonList))
# 
# for(j in 1:length(InputList[[i]]$TaxonList)){
#   TaxVect[InputList[[i]]$TaxonList[[j]]] <- j
# }
# 
# WorkStruct <- InputList[[i]]$PrinCurveStruct$ListProc[[length(InputList[[i]]$PrinCurveStruct$ListProc)]]
# SampleReord <- lapply(as.list(1:100), function(i){sample(TaxVect)})
# PvVect <- rep(NA, nrow(InputList[[i]]$Expression))
# 
# DoStuff <- function(gId) {
#   Exp <- unlist(InputList[[i]]$Expression[gId,names(WorkStruct$CellsPT)])
#   Base <- median(unlist(lapply(split(Exp, f = factor(TaxVect)), mad)))
#   
#   MedVect <- sapply(SampleReord, function(x){
#     median(unlist(lapply(split(Exp, f = factor(x)), mad)))
#   })
#   return(wilcox.test(MedVect - Base, alternative = "greater")$p.value)
# }
# 
# cl <- parallel::makeCluster(4, type = "FORK")
# parallel::clusterExport(cl, varlist = c("InputList", "WorkStruct", "TaxVect", "SampleReord", "i"), envir = environment())
# 
# AllPV.Sasa <- parallel::parLapply(cl, as.list(1:nrow(InputList[[i]]$Expression)), DoStuff)
# 
# 
# AllPV.Sasa <- unlist(AllPV.Sasa)
# names(AllPV.Sasa) <- rownames(InputList[[i]]$Expression)
# 
# 
# write_rds(AllPV.Sasa, path = "AllPV.Sasa.rds")
# 









# 
# INCCGO <- (names(AllPV) %in% MouseGenes_GOCellCycle)
# boxplot(AllPV ~ INCCGO, log='y')
# 
# 

























ProjectOnCircle <- function(Points, Edges, Nodes, Categories, Title, ExpValues, PCACenter = TRUE) {
  
  if(PCACenter){
    ScaledNodes <- scale(Nodes, center = PCACenter, scale = FALSE)
    Centers <- attr(ScaledNodes, "scaled:center")
    Nodes <- ScaledNodes
  }
  
  PCAPrGraph <- prcomp(Nodes, retx = TRUE, center = FALSE, scale. = FALSE)
  VarExp <- PCAPrGraph$sdev[1:2]^2/sum(PCAPrGraph$sdev^2)
  
  if(is.null(Categories)){
    Categories <- rep("NoG", nrow(Points))
  }
  
  if(PCACenter){
    Points <- scale(Points, center = Centers, scale = FALSE)
  }
  
  RotatedPoints <- Points %*% PCAPrGraph$rotation[,1:2]
  
  RotatedData <- cbind(RotatedPoints,
                       as.character(Categories),
                       ExpValues[match(rownames(RotatedPoints), names(ExpValues))])
  
  
  colnames(RotatedData) <- c("PC1", "PC2", "Cat", "Exp")
  
  RotatedData.DF <- data.frame(RotatedData)
  RotatedData.DF$PC1 <- as.numeric(as.character(RotatedData.DF$PC1))
  RotatedData.DF$PC2 <- as.numeric(as.character(RotatedData.DF$PC2))
  RotatedData.DF$Exp <- as.numeric(as.character(RotatedData.DF$Exp))
  
  p <- ggplot2::ggplot(data.frame(RotatedData.DF), ggplot2::aes(x = PC1, y = PC2, colour = Exp, shape = Cat)) +
    ggplot2::scale_color_continuous(low = "blue", high = "red")
  
  p <- p + ggplot2::geom_point(data = data.frame(PCAPrGraph$x[,1:2]), mapping = ggplot2::aes(x=PC1, y=PC2), inherit.aes = FALSE) +
    ggplot2::labs(title = Title, x = paste("PC1 -", signif(100*VarExp[1], 4), "%"), y = paste("PC2 -", signif(100*VarExp[2], 4), "%"))
  
  for(j in 1:nrow(Edges)){
    p <- p + ggplot2::geom_path(data = data.frame(PCAPrGraph$x[Edges[j,],1:2]),
                                mapping = ggplot2::aes(x = PC1, y = PC2), inherit.aes = FALSE)
  }
  
  p <- p + ggplot2::geom_point()
  
  return(p)
  
}















# Genes.Buett <- sort(unique(c(rownames(FinalStructure$NodesExp), rownames(FinalStructure$NodesExp))))

# gName <- sample(MouseGenes_GOCellCycle, 1)
# 
# if(gName %in% rownames(FinalStructure.Buett$NodesExp)){
#   p.Buett <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Buett$NodesPT,
#                                                y=FinalStructure.Buett$NodesExp[gName,]),
#                              mapping = ggplot2::aes(x = x, y = y, color="PC")) +
#     ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(gName, " - Buettner et al.", sep ="")) +
#     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
#     ggplot2::geom_rect(data = FinalStructure.Buett$RecCoord,
#                        mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
#                        ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
#     ggplot2::geom_point(data = data.frame(x=FinalStructure.Buett$CellsPT[ReOrd.Buett],
#                                           y=FinalStructure.Buett$CellExp[gName,]),
#                         mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
#     ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))
#   
# } else {
#   p.Buett <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Buett$NodesPT,
#                                                y=FinalStructure.Buett$NodesExp[1,]),
#                              mapping = ggplot2::aes(x = x, y = y, color="PC"))
# }
# 
# if(gName %in% rownames(FullExpData.Buett)){
#   p.Buett.2 <- ProjectOnCircle(
#     Points = ExpDataset.Buett$ListInfo[[length(ExpDataset.Buett$ListInfo)]]$FinalStruct$Data,
#     Edges = ExpDataset.Buett$ListInfo[[length(ExpDataset.Buett$ListInfo)]]$FinalStruct$IntGrahs[[
#       length(ExpDataset.Buett$ListInfo[[length(ExpDataset.Buett$ListInfo)]]$FinalStruct$IntGrahs)
#       ]]$Edges,
#     Nodes = ExpDataset.Buett$ListInfo[[length(ExpDataset.Buett$ListInfo)]]$FinalStruct$IntGrahs[[
#       length(ExpDataset.Buett$ListInfo[[length(ExpDataset.Buett$ListInfo)]]$FinalStruct$IntGrahs)
#       ]]$Nodes,
#     Categories = ExpDataset.Buett$ListInfo[[length(ExpDataset.Buett$ListInfo)]]$FinalStruct$Categories,
#     Title = gName,
#     ExpValues = unlist(FullExpData.Buett[gName,])
#   )
# } else {
#   p.Buett.2 <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Buett$NodesPT,
#                                                y=FinalStructure.Buett$NodesExp[1,]),
#                              mapping = ggplot2::aes(x = x, y = y, color="PC"))
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
# if(gName %in% rownames(FinalStructure.Sasa$NodesExp)){
#   p.Sasa <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Sasa$NodesPT,
#                                               y=FinalStructure.Sasa$NodesExp[gName,]),
#                             mapping = ggplot2::aes(x = x, y = y, color="PC")) +
#     ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(gName, " - Sasagawa et al.", sep ="")) +
#     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
#     ggplot2::geom_rect(data = FinalStructure.Sasa$RecCoord,
#                        mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
#                        ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
#     ggplot2::geom_point(data = data.frame(x=FinalStructure.Sasa$CellsPT[ReOrd.Sasa],
#                                           y=FinalStructure.Sasa$CellExp[gName,]),
#                         mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
#     ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))
# } else {
#   p.Sasa <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Sasa$NodesPT,
#                                               y=FinalStructure.Sasa$NodesExp[1,]),
#                             mapping = ggplot2::aes(x = x, y = y, color="PC"))
# }
#   
# if(gName %in% rownames(FullExpData.Sasa)){
#   p.Sasa.2 <- ProjectOnCircle(
#     Points = ExpDataset.Sasa$ListInfo[[length(ExpDataset.Sasa$ListInfo)]]$FinalStruct$Data,
#     Edges = ExpDataset.Sasa$ListInfo[[length(ExpDataset.Sasa$ListInfo)]]$FinalStruct$IntGrahs[[
#       length(ExpDataset.Sasa$ListInfo[[length(ExpDataset.Sasa$ListInfo)]]$FinalStruct$IntGrahs)
#       ]]$Edges,
#     Nodes = ExpDataset.Sasa$ListInfo[[length(ExpDataset.Sasa$ListInfo)]]$FinalStruct$IntGrahs[[
#       length(ExpDataset.Sasa$ListInfo[[length(ExpDataset.Sasa$ListInfo)]]$FinalStruct$IntGrahs)
#       ]]$Nodes,
#     Categories = ExpDataset.Sasa$ListInfo[[length(ExpDataset.Sasa$ListInfo)]]$FinalStruct$Categories,
#     Title = gName,
#     ExpValues = unlist(FullExpData.Sasa[gName,])
#   )
# } else {
#   p.Sasa2 <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Sasa$NodesPT,
#                                               y=FinalStructure.Sasa$NodesExp[1,]),
#                             mapping = ggplot2::aes(x = x, y = y, color="PC"))
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
# if(gName %in% rownames(FinalStructure.Kowa$NodesExp)){
#   p.Kowa <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Kowa$NodesPT,
#                                               y=FinalStructure.Kowa$NodesExp[gName,]),
#                             mapping = ggplot2::aes(x = x, y = y, color="PC")) +
#     ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(gName, " - Kowalczyk et al.", sep ="")) +
#     ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
#     ggplot2::geom_rect(data = FinalStructure.Kowa$RecCoord,
#                        mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
#                        ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
#     ggplot2::geom_point(data = data.frame(x=FinalStructure.Kowa$CellsPT[ReOrd.Kowa],
#                                           y=FinalStructure.Kowa$CellExp[gName,]),
#                         mapping = ggplot2::aes(x=x, y=y, color="Data"), inherit.aes = FALSE, alpha=.5) +
#     ggplot2::geom_point() + ggplot2::geom_line() + ggplot2::scale_color_manual(values = c("blue", "black"))
# } else {
#   p.Kowa <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Kowa$NodesPT,
#                                               y=FinalStructure.Kowa$NodesExp[1,]),
#                             mapping = ggplot2::aes(x = x, y = y, color="PC"))
# }
# 
# if(gName %in% rownames(FullExpData.Sasa)){
#   p.Kowa.2 <- ProjectOnCircle(
#     Points = ExpDataset.Kowa$ListInfo[[length(ExpDataset.Kowa$ListInfo)]]$FinalStruct$Data,
#     Edges = ExpDataset.Kowa$ListInfo[[length(ExpDataset.Kowa$ListInfo)]]$FinalStruct$IntGrahs[[
#       length(ExpDataset.Kowa$ListInfo[[length(ExpDataset.Kowa$ListInfo)]]$FinalStruct$IntGrahs)
#       ]]$Edges,
#     Nodes = ExpDataset.Kowa$ListInfo[[length(ExpDataset.Kowa$ListInfo)]]$FinalStruct$IntGrahs[[
#       length(ExpDataset.Kowa$ListInfo[[length(ExpDataset.Kowa$ListInfo)]]$FinalStruct$IntGrahs)
#       ]]$Nodes,
#     Categories = ExpDataset.Kowa$ListInfo[[length(ExpDataset.Kowa$ListInfo)]]$FinalStruct$Categories,
#     Title = gName,
#     ExpValues = unlist(FullExpData.Kowa[gName,])
#   )
# } else {
#   p.Kowa.2 <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Kowa$NodesPT,
#                                               y=FinalStructure.Kowa$NodesExp[1,]),
#                             mapping = ggplot2::aes(x = x, y = y, color="PC"))
# }
# 
# 
# 
# 
# 
# 
# 
# 













# gridExtra::grid.arrange(p.Buett, p.Kowa, p.Sasa, ncol=2)
# 
# gridExtra::grid.arrange(p.Buett, p.Buett.2, ncol=1)
# gridExtra::grid.arrange(p.Kowa, p.Kowa.2, ncol=1)
# gridExtra::grid.arrange(p.Sasa, p.Sasa.2, ncol=1)

# Sys.sleep(1)
# 
# gridExtra::grid.arrange(p.Buett, p.Buett.2,
#                         p.Kowa, p.Kowa.2,
#                         p.Sasa, p.Sasa.2, ncol=2)









# print(p.Buett)
# print(p.Sasa)
# print(p.Kowa)




SumPeack <- sapply(as.list(rownames(FinalStructure.Kowa$NodesExp)), function(gName){
  PosIdxs <- (1:ncol(FinalStructure.Kowa$NodesExp))[FinalStructure.Kowa$NodesExp[gName,]
                                                    >= quantile(FinalStructure.Kowa$NodesExp[gName,], .9)]

  if(length(PosIdxs) == 0){
    return(rep(FALSE, length(FinalStructure.Kowa$StageOnNodes)+1))
  } else {
    c(sapply(FinalStructure.Kowa$StageOnNodes, function(ids){
      any(ids %in% which.max(FinalStructure.Kowa$NodesExp[gName,]))
    }),all(min(PosIdxs):max(PosIdxs) %in% PosIdxs))
  }
  
})

rownames(SumPeack)[nrow(SumPeack)] <- "Strong"
colnames(SumPeack) <- rownames(FinalStructure.Kowa$NodesExp)

dim(SumPeack)

SumPeack <- SumPeack[,SumPeack["Strong", ]]

dim(SumPeack)

GenesPhase.Kowa <- apply(SumPeack[-nrow(SumPeack),], 1, which)

barplot(unlist(lapply(GenesPhase.Kowa, length)), las=2)










SumPeack <- sapply(as.list(rownames(FinalStructure.Sasa$NodesExp)), function(gName){
  PosIdxs <- (1:ncol(FinalStructure.Sasa$NodesExp))[FinalStructure.Sasa$NodesExp[gName,]
                                                    >= quantile(FinalStructure.Sasa$NodesExp[gName,], .9)]
  
  if(length(PosIdxs) == 0){
    return(rep(FALSE, length(FinalStructure.Sasa$StageOnNodes)+1))
  } else {
    c(sapply(FinalStructure.Sasa$StageOnNodes, function(ids){
      any(ids %in% which.max(FinalStructure.Sasa$NodesExp[gName,]))
    }),all(min(PosIdxs):max(PosIdxs) %in% PosIdxs))
  }

})

rownames(SumPeack)[nrow(SumPeack)] <- "Strong"
colnames(SumPeack) <- rownames(FinalStructure.Sasa$NodesExp)

dim(SumPeack)

SumPeack <- SumPeack[,SumPeack["Strong", ]]

dim(SumPeack)

GenesPhase.Sasa <- apply(SumPeack[-nrow(SumPeack),], 1, which)

barplot(unlist(lapply(GenesPhase.Sasa, length)), las=2)








SumPeack <- sapply(as.list(rownames(FinalStructure.Buett$NodesExp)), function(gName){
  PosIdxs <- (1:ncol(FinalStructure.Buett$NodesExp))[FinalStructure.Buett$NodesExp[gName,]
                                                    >= quantile(FinalStructure.Buett$NodesExp[gName,], .9)]
  
  if(length(PosIdxs) == 0){
    return(rep(FALSE, length(FinalStructure.Buett$StageOnNodes)+1))
  } else {
    c(sapply(FinalStructure.Buett$StageOnNodes, function(ids){
      any(ids %in% which.max(FinalStructure.Buett$NodesExp[gName,]))
    }),all(min(PosIdxs):max(PosIdxs) %in% PosIdxs))
  }
})

rownames(SumPeack)[nrow(SumPeack)] <- "Strong"
colnames(SumPeack) <- rownames(FinalStructure.Buett$NodesExp)

dim(SumPeack)

SumPeack <- SumPeack[,SumPeack["Strong", ]]

dim(SumPeack)

GenesPhase.Buett <- apply(SumPeack[-nrow(SumPeack),], 1, which)



barplot(unlist(lapply(GenesPhase.Buett, length)), las=2)








intersect(c(names(GenesPhase.Buett$G1), names(GenesPhase.Buett$`G1+S`)),
          c(names(GenesPhase.Sasa$S_G1), names(GenesPhase.Sasa$`S_G1+S_S`)))

intersect(c(names(GenesPhase.Buett$S), names(GenesPhase.Buett$`G1+S`), names(GenesPhase.Buett$`S+G2M`)),
          c(names(GenesPhase.Sasa$S_S), names(GenesPhase.Sasa$`S_G1+S_S`), names(GenesPhase.Sasa$`S_S+S_G2/M`)))

intersect(c(names(GenesPhase.Buett$G2M), names(GenesPhase.Buett$`G2M+G1`), names(GenesPhase.Buett$`S+G2M`)),
          c(names(GenesPhase.Sasa$`S_G2/M`), names(GenesPhase.Sasa$`S_G2/M+S_G1`), names(GenesPhase.Sasa$`S_S+S_G2/M`)))







G1Shared <- intersect(c(names(GenesPhase.Buett$G1), names(GenesPhase.Buett$`G1+S`)),
                      c(names(GenesPhase.Kowa$`G0+G1(early)`), names(GenesPhase.Kowa$`G1(early)+G1(late)`),
                        names(GenesPhase.Kowa$`G1(late)`), names(GenesPhase.Kowa$`G1(late)+S`)))

SShared <- intersect(c(names(GenesPhase.Buett$S), names(GenesPhase.Buett$`G1+S`), names(GenesPhase.Buett$`S+G2M`)),
                     c(names(GenesPhase.Kowa$`G1(late)+S`), names(GenesPhase.Kowa$S), names(GenesPhase.Kowa$`S+G2/M`)))

G2Shared <- intersect(c(names(GenesPhase.Buett$G2M), names(GenesPhase.Buett$`G2M+G1`), names(GenesPhase.Buett$`S+G2M`)),
                      c(names(GenesPhase.Kowa$`S+G2/M`), names(GenesPhase.Kowa$`G2/M`), names(GenesPhase.Kowa$`G2/M+G0`)))






G1Any <- union(c(names(GenesPhase.Buett$G1), names(GenesPhase.Buett$`G1+S`)),
               c(names(GenesPhase.Kowa$`G0+G1(early)`), names(GenesPhase.Kowa$`G1(early)+G1(late)`),
                 names(GenesPhase.Kowa$`G1(late)`), names(GenesPhase.Kowa$`G1(late)+S`)))

SAny <- union(c(names(GenesPhase.Buett$S), names(GenesPhase.Buett$`G1+S`), names(GenesPhase.Buett$`S+G2M`)),
              c(names(GenesPhase.Kowa$`G1(late)+S`), names(GenesPhase.Kowa$S), names(GenesPhase.Kowa$`S+G2/M`)))

G2Any <- union(c(names(GenesPhase.Buett$G2M), names(GenesPhase.Buett$`G2M+G1`), names(GenesPhase.Buett$`S+G2M`)),
               c(names(GenesPhase.Kowa$`S+G2/M`), names(GenesPhase.Kowa$`G2/M`), names(GenesPhase.Kowa$`G2/M+G0`)))












GMTSel <- rRoma::SelectFromMSIGdb("WHITFIELD_CELL_CYCLE")













gName <- sample(MouseGenes_GOCellCycle, 1)


AllCells <- rbind(data.frame(x=FinalStructure.Buett$CellsPT[ReOrd.Buett]/max(FinalStructure.Buett$CellsPT),
                             y=unlist(FullExpData.Buett[gName, names(FinalStructure.Buett$CellsPT[ReOrd.Buett])]),
                             Org=rep("Buet", length(ReOrd.Buett))),
                  data.frame(x=FinalStructure.Sasa$CellsPT[ReOrd.Sasa]/max(FinalStructure.Sasa$CellsPT),
                             y=unlist(FullExpData.Sasa[gName, names(FinalStructure.Sasa$CellsPT[ReOrd.Sasa])]),
                             Org=rep("Sasa", length(ReOrd.Sasa))),
                  data.frame(x=FinalStructure.Kowa$CellsPT[ReOrd.Kowa]/max(FinalStructure.Kowa$CellsPT),
                             y=unlist(FullExpData.Kowa[gName, names(FinalStructure.Kowa$CellsPT[ReOrd.Kowa])]),
                             Org=rep("Kowa", length(ReOrd.Kowa)))
                  )


ReNorm.Buett <- FinalStructure.Buett$RecCoord
ReNorm.Buett[,1:3] <- ReNorm.Buett[,1:3]/max(FinalStructure.Buett$CellsPT)

ReNorm.Sasa <- FinalStructure.Sasa$RecCoord
ReNorm.Sasa[,1:3] <- ReNorm.Sasa[,1:3]/max(FinalStructure.Sasa$CellsPT)

ReNorm.Kowa <- FinalStructure.Kowa$RecCoord
ReNorm.Kowa[,1:3] <- ReNorm.Kowa[,1:3]/max(FinalStructure.Kowa$CellsPT)

AllRect <- rbind(cbind(ReNorm.Buett, Org = rep("Buet", nrow(FinalStructure.Buett$RecCoord))),
                 cbind(ReNorm.Sasa, Org = rep("Sasa", nrow(FinalStructure.Sasa$RecCoord))),
                 cbind(ReNorm.Kowa, Org = rep("Kowa", nrow(FinalStructure.Kowa$RecCoord))))

AllRect$Stage <- as.character(AllRect$Stage)



AllRect$Stage[AllRect$Stage == "G2/M"] <- "G2M"
AllRect$Stage[AllRect$Stage == "G2/M+G0"] <- "G2M+G0"
AllRect$Stage[AllRect$Stage == "S_G1"] <- "G1"
AllRect$Stage[AllRect$Stage == "S_S"] <- "G1"
AllRect$Stage[AllRect$Stage == "S_G2/M"] <- "G2M"


ggplot2::ggplot(data = AllCells, mapping = ggplot2::aes(x=x, y=y)) +
  ggplot2::labs(x = "Pseudotime", y="Gene expression", title = gName) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::geom_rect(data = AllRect, mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                     ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) +
  ggplot2::geom_smooth(span = 0.5) + ggplot2::geom_point(alpha=.5) +
  ggplot2::facet_grid(facets = Org ~ ., scales = "free") 
  







p.Buett <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Buett$CellsPT[ReOrd.Buett],
                                             y=unlist(FullExpData.Buett[gName, names(FinalStructure.Buett$CellsPT[ReOrd.Buett])])),
                           mapping = ggplot2::aes(x=x, y=y, color="Data")) +
  ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(gName, " - Buettner et al.", sep ="")) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::geom_rect(data = FinalStructure.Buett$RecCoord,
                     mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                     ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
  ggplot2::geom_point(alpha=.5) + ggplot2::geom_smooth(span = 0.5) + ggplot2::scale_color_manual(values = c("blue", "black"))
# print(p.Buett)


p.Sasa <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Sasa$CellsPT[ReOrd.Sasa],
                                             y=unlist(FullExpData.Sasa[gName, names(FinalStructure.Sasa$CellsPT[ReOrd.Sasa])])),
                           mapping = ggplot2::aes(x=x, y=y, color="Data")) +
  ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(gName, " - Sasagawa et al.", sep ="")) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::geom_rect(data = FinalStructure.Sasa$RecCoord,
                     mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                     ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
  ggplot2::geom_point(alpha=.5) + ggplot2::geom_smooth(span = 0.5) + ggplot2::scale_color_manual(values = c("blue", "black"))
# print(p.Sasa)




p.Kowa <- ggplot2::ggplot(data = data.frame(x=FinalStructure.Kowa$CellsPT[ReOrd.Kowa],
                                            y=unlist(FullExpData.Kowa[gName, names(FinalStructure.Kowa$CellsPT[ReOrd.Kowa])])),
                          mapping = ggplot2::aes(x=x, y=y, color="Data")) +
  ggplot2::labs(x = "Pseudotime", y="Gene expression", title = paste(gName, " - Kowalczyk et al.", sep ="")) +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
  ggplot2::geom_rect(data = FinalStructure.Kowa$RecCoord,
                     mapping = ggplot2::aes(fill=Stage, xmin=Min, xmax=Max),
                     ymin = -Inf, ymax = Inf, inherit.aes = FALSE, alpha=.4) + 
  ggplot2::geom_point(alpha=.5) + ggplot2::geom_smooth(span = 0.5) + ggplot2::scale_color_manual(values = c("blue", "black"))
# print(p.Kowa)





gridExtra::grid.arrange(p.Buett, p.Kowa, p.Sasa, ncol=2)

# gridExtra::grid.arrange(p.Buett, p.Sasa, ncol=1)






















library(VennDiagram)
grid.newpage()


draw.pairwise.venn(area1 = length(Genes.Buett), area2 = length(MouseGenes_GOCellCycle),
                   cross.area = length(intersect(Genes.Buett, MouseGenes_GOCellCycle)),
                   fill = c("red", "blue"), cex = 3,
                   category = c("Buettner et al.", "GO"))
grid.newpage()



draw.pairwise.venn(area1 = length(Genes.Sasa), area2 = length(MouseGenes_GOCellCycle),
                   cross.area = length(intersect(Genes.Sasa, MouseGenes_GOCellCycle)),
                   fill = c("red", "blue"), cex = 3,
                   category = c("Sasagawa et al.", "GO"))
grid.newpage()


draw.pairwise.venn(area1 = length(Genes.Kowa), area2 = length(MouseGenes_GOCellCycle),
                   cross.area = length(intersect(Genes.Kowa, MouseGenes_GOCellCycle)),
                   fill = c("red", "blue"), cex = 3,
                   category = c("Kowalczyk et al.", "GO"))
grid.newpage()



draw.pairwise.venn(area1 = length(unique(c(Genes.Kowa, Genes.Buett, Genes.Sasa))),
                   area2 = length(MouseGenes_GOCellCycle),
                   cross.area = length(intersect(unique(c(Genes.Kowa, Genes.Buett, Genes.Sasa)), MouseGenes_GOCellCycle)),
                   fill = c("red", "blue"), cex = 3,
                   category = c("B+S+K", "GO"))
grid.newpage()






draw.pairwise.venn(area1 = length(Genes.Buett), area2 = length(Genes.Sasa),
                   cross.area = length(intersect(Genes.Buett, Genes.Sasa)),
                   fill = c("red", "blue"), cex = 3,
                   category = c("Buettner et al.", "Sasagawa et al."))
grid.newpage()



draw.triple.venn(area1 = length(Genes.Buett), area2 = length(MouseGenes_GOCellCycle), area3 = length(Genes.Sasa),
                 n12 = length(intersect(Genes.Buett, MouseGenes_GOCellCycle)),
                 n23 = length(intersect(Genes.Sasa, MouseGenes_GOCellCycle)),
                 n13 = length(intersect(Genes.Buett, Genes.Sasa)),
                 n123 = length(intersect(Genes.Buett, intersect(MouseGenes_GOCellCycle, Genes.Sasa))),
                 category = c("Buettner et al.", "GO", "Sasagawa et al."),
                 fill = c("red", "white", "blue"), cex = 3)
grid.newpage()





draw.triple.venn(area1 = length(Genes.Buett), area2 = length(Genes.Kowa), area3 = length(Genes.Sasa),
                 n12 = length(intersect(Genes.Buett, Genes.Kowa)),
                 n23 = length(intersect(Genes.Sasa, Genes.Kowa)),
                 n13 = length(intersect(Genes.Buett, Genes.Sasa)),
                 n123 = length(intersect(Genes.Buett, intersect(Genes.Kowa, Genes.Sasa))),
                 category = c("Buettner et al.", "Kowalczyk et al.", "Sasagawa et al."),
                 fill = c("red", "green", "blue"), cex = 3)
grid.newpage()



draw.quad.venn(area1 = length(Genes.Buett), area2 = length(Genes.Kowa),
               area3 = length(Genes.Sasa), area4 = length(MouseGenes_GOCellCycle),
               n12 = length(intersect(Genes.Buett, Genes.Kowa)),
               n23 = length(intersect(Genes.Sasa, Genes.Kowa)),
               n13 = length(intersect(Genes.Buett, Genes.Sasa)),
               n14 = length(intersect(Genes.Buett, MouseGenes_GOCellCycle)),
               n24 = length(intersect(MouseGenes_GOCellCycle, Genes.Kowa)),
               n34 = length(intersect(Genes.Sasa, MouseGenes_GOCellCycle)),
               n123 = length(intersect(Genes.Buett, intersect(Genes.Kowa, Genes.Sasa))),
               n124 = length(intersect(Genes.Buett, intersect(Genes.Kowa, MouseGenes_GOCellCycle))),
               n134 = length(intersect(Genes.Buett, intersect(Genes.Sasa, MouseGenes_GOCellCycle))),
               n234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, MouseGenes_GOCellCycle))),
               n1234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, intersect(Genes.Buett, MouseGenes_GOCellCycle)))),
               category = c("Buettner et al.", "Kowalczyk et al.", "Sasagawa et al.", "GO"),
               fill = c("red", "green", "blue", "white"), cex = 3)
grid.newpage()
































X <- FinalStructure.Buett$CellsPT[ReOrd.Buett]/max(FinalStructure.Buett$NodesPT)
Xfit <- c(X, X + 1, X - 1)

cl <- parallel::makeCluster(3, type = "FORK")

parallel::clusterExport(cl=cl, varlist=c("Xfit"), envir = environment())

tictoc::tic()
AllModelFit.Buett <- parallel::parApply(cl, FullExpData.Buett[,names(FinalStructure.Buett$CellsPT[ReOrd.Buett])], 1, function(y){
  
  Yfit <- rep(unlist(y), 3)
  Loess <- loess(Yfit ~ Xfit, span = .5)
  predict(Loess, newdata = data.frame(Xfit = unique(X)))
  
})
tictoc::toc()

parallel::stopCluster(cl)

DiffMat <- t(AllModelFit.Buett) - FullExpData.Buett[,names(FinalStructure.Buett$CellsPT[ReOrd.Buett])]

BaseRes <- apply(FullExpData.Buett[,names(FinalStructure.Buett$CellsPT[ReOrd.Buett])], 1, function(x){
  sum((x - mean(x))^2)
})

ResRat.Buett <- BaseRes/rowSums(DiffMat^2)

SumPeack <- sapply(as.list(colnames(AllModelFit.Buett)), function(gName){
  PosIdxs <- (1:nrow(AllModelFit.Buett))[AllModelFit.Buett[,gName]
                                                     >= quantile(AllModelFit.Buett[,gName], .9)]
  
  # Need to associate cells and stages
  
  if(length(PosIdxs) == 0){
    return(rep(FALSE, length(FinalStructure.Buett$StageOnNodes)+1))
  } else {
    c(sapply(FinalStructure.Buett$StageOnNodes, function(ids){
      any(ids %in% which.max(FinalStructure.Buett$NodesExp[gName,]))
    }),all(min(PosIdxs):max(PosIdxs) %in% PosIdxs))
  }
})

rownames(SumPeack)[nrow(SumPeack)] <- "Strong"
colnames(SumPeack) <- rownames(FinalStructure.Buett$NodesExp)

dim(SumPeack)

SumPeack <- SumPeack[,SumPeack["Strong", ]]

dim(SumPeack)

GenesPhase.Buett <- apply(SumPeack[-nrow(SumPeack),], 1, which)
















X <- FinalStructure.Kowa$CellsPT[ReOrd.Kowa]/max(FinalStructure.Kowa$NodesPT)
Xfit <- c(X, X + 1, X - 1)

cl <- parallel::makeCluster(3, type = "FORK")

parallel::clusterExport(cl=cl, varlist=c("Xfit"), envir = environment())

tictoc::tic()
AllModelFit.Kowa <- parallel::parApply(cl, FullExpData.Kowa[,names(FinalStructure.Kowa$CellsPT[ReOrd.Kowa])], 1, function(y){
  
  Yfit <- rep(unlist(y), 3)
  Loess <- loess(Yfit ~ Xfit, span = .5)
  predict(Loess, newdata = data.frame(Xfit = X))
  
})
tictoc::toc()

parallel::stopCluster(cl)

DiffMat <- t(AllModelFit.Kowa) - FullExpData.Kowa[,names(FinalStructure.Kowa$CellsPT[ReOrd.Kowa])]


BaseRes <- apply(FullExpData.Kowa[,names(FinalStructure.Kowa$CellsPT[ReOrd.Kowa])], 1, function(x){
  sum((x - mean(x))^2)
})


ResRat.Kowa <- BaseRes/rowSums(DiffMat^2)












X <- FinalStructure.Sasa$CellsPT[ReOrd.Sasa]/max(FinalStructure.Sasa$NodesPT)
Xfit <- c(X, X + 1, X - 1)

cl <- parallel::makeCluster(3, type = "FORK")

parallel::clusterExport(cl=cl, varlist=c("Xfit"), envir = environment())

tictoc::tic()
AllModelFit.Sasa <- parallel::parApply(cl, FullExpData.Sasa[,names(FinalStructure.Sasa$CellsPT[ReOrd.Sasa])], 1, function(y){
  
  Yfit <- rep(unlist(y), 3)
  Loess <- loess(Yfit ~ Xfit, span = .5)
  predict(Loess, newdata = data.frame(Xfit = X))
  
})
tictoc::toc()

parallel::stopCluster(cl)

DiffMat <- t(AllModelFit.Sasa) - FullExpData.Sasa[,names(FinalStructure.Sasa$CellsPT[ReOrd.Sasa])]


BaseRes <- apply(FullExpData.Sasa[,names(FinalStructure.Sasa$CellsPT[ReOrd.Sasa])], 1, function(x){
  sum((x - mean(x))^2)
})


ResRat.Sasa <- BaseRes/rowSums(DiffMat^2)




AllGenes <- unique(c(names(ResRat.Buett), names(ResRat.Kowa), names(ResRat.Sasa)))





AcrossDS <- rbind(ResRat.Sasa[AllGenes],
                  ResRat.Kowa[AllGenes],
                  ResRat.Buett[AllGenes])


AcrossDS[,which(apply(AcrossDS < 1, 2, all, na.rm=TRUE))[1:50]]



AllNames <- names(which(apply(AcrossDS < 1, 2, all, na.rm=TRUE)))






AllGoodGenes <- names(which(apply(AcrossDS[, colSums(is.na(AcrossDS))==0] < .98, 2, any)))


draw.quad.venn(area1 = length(Genes.Buett), area2 = length(Genes.Kowa),
               area3 = length(Genes.Sasa), area4 = length(AllGoodGenes),
               n12 = length(intersect(Genes.Buett, Genes.Kowa)),
               n23 = length(intersect(Genes.Sasa, Genes.Kowa)),
               n13 = length(intersect(Genes.Buett, Genes.Sasa)),
               n14 = length(intersect(Genes.Buett, AllGoodGenes)),
               n24 = length(intersect(AllGoodGenes, Genes.Kowa)),
               n34 = length(intersect(Genes.Sasa, AllGoodGenes)),
               n123 = length(intersect(Genes.Buett, intersect(Genes.Kowa, Genes.Sasa))),
               n124 = length(intersect(Genes.Buett, intersect(Genes.Kowa, AllGoodGenes))),
               n134 = length(intersect(Genes.Buett, intersect(Genes.Sasa, AllGoodGenes))),
               n234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, AllGoodGenes))),
               n1234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, intersect(Genes.Buett, AllGoodGenes)))),
               category = c("Buettner et al.", "Kowalczyk et al.", "Sasagawa et al.", "Smooth"),
               fill = c("red", "green", "blue", "white"), cex = 3)
grid.newpage()





"Ccnd2"   "Cdkn2c"   
   "Npat"    "Cdc6"    "Tex14"   "Pkp4"   "Ccne1"   "Fam175a" "Egf"        "Mcm3"   


