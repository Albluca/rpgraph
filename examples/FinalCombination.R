library(readr)
library(Pricecycle)

# Load data

Data.Sasa <- read_rds("~/Google Drive/Datasets/Sasagawa et al - Murine Stem Cells (EB5 cell line)/rPG_Last.rds")
Data.Kowa <- read_rds("~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.rds")
Data.Buet <- read_rds("~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/rPG_Last.rds")

# Plot data and construct auxiliary structures

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
                              SelThr = .3,
                              ComputeOverlaps = TRUE,
                              ExpData = TargetStruct$FiltExp,
                              RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                              PCACenter = TargetStruct$PCAData$center,
                              PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                              OrderOnCat = TRUE,
                              SmoothPoints = 0, MinCellPerNode = 2)



# Create structure to use thge dashboard


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
list(Name = "Sasagawa et al", Expression = Data.Sasa$ExpMat,
     Categories = Data.Sasa$Cats,
     OrderedData = Proc.Exp.Sasa,
     PGStruct = Data.Sasa$Analysis$PGStructs[[length(Data.Sasa$Analysis$PGStructs)]],
     TaxonList = Data.Sasa$Analysis$PGStructs[[length(Data.Sasa$Analysis$PGStructs)]]$TaxonList
)
)


# Launch dashboard

CompareAcrossData(InputList, NULL)

# CatOrder = c("G0", "G0+G1", "G1", "G1+S", "S", "S+G2M", "G2M", "G2M+G0", "G1M+G1")






# Proc.Exp.Buet$StageOnNodes
#
#
# WorkPG <- Data.Sasa$Analysis$PGStructs[[length(Data.Sasa$Analysis$PGStructs)]]
# TaxList <- WorkPG$TaxonList[[length(WorkPG$TaxonList)]]
# TaxVect <- rep(NA, nrow(WorkPG$Data))
# names(TaxVect) <- rownames(WorkPG$Data)
#
# lapply(1:length(TaxList), function(j){
#   TaxVect[ TaxList[[j]] ] <<- j
# })
#
#
# GeneExprMat.DF <- data.frame(t(Data.Sasa$ExpMat[, names(Proc.Exp.Sasa$CellsPT)]))
# colnames(GeneExprMat.DF) <- rownames(Data.Sasa$ExpMat)
#
# DoStuff <- function(TaxVect){
#
#   SelCell <- names(TaxVect)
#
#   GeneExprMat.DF.Split <- split(GeneExprMat.DF[SelCell, ], TaxVect[SelCell])
#
#   GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)
#
#   GeneExprMat.DF.MeanRemoved <-
#     lapply(as.list(1:length(GeneExprMat.DF.Split)), function(i){
#       apply(GeneExprMat.DF.Split[[i]], 2, sd)/GeneExprMat.DF.Split.Mean[[i]]
#     })
#
#   GeneExprMat.DF.MeanRemoved.All <- NULL
#   for(i in 1:length(GeneExprMat.DF.MeanRemoved)){
#     GeneExprMat.DF.MeanRemoved.All <- rbind(GeneExprMat.DF.MeanRemoved.All, GeneExprMat.DF.MeanRemoved[[i]])
#   }
#
#   return(apply(abs(GeneExprMat.DF.MeanRemoved.All), 2, mean, na.rm = TRUE))
#
# }
#
# Base <- DoStuff(TaxVect)
# Base <- Base[is.finite(Base)]
#
#
#
#
# plot(Base, apply(GeneExprMat.DF[,names(Base)], 2, mean), xlab = "Mean CV", ylab = "Mean Expression")
# points(Base[Genes.Sasa], apply(GeneExprMat.DF[,Genes.Sasa], 2, mean), col='red')
#
#
# BinMat <- t(GeneExprMat.DF) > apply(GeneExprMat.DF[,], 2, quantile, probs = .9)
#
# pheatmap::pheatmap(1*BinMat[1:100, ])
#
# pheatmap::pheatmap(1*BinMat[1:2, ])
#
#
# plot(apply(GeneExprMat.DF[,names(Base)], 2, sd), apply(GeneExprMat.DF[,names(Base)], 2, mean), xlab = "Mean CV", ylab = "Mean Expression")


DataStruct = Data.Buet
ProcStruct = Proc.Exp.Buet
FiltMax = .5
Thr = .7

GetGenesWithPeaks <- function(DataStruct, ProcStruct, FiltMax = .2, Thr = .9) {

  NormExp <- ProcStruct$NodesExp
  NormExp[NormExp < FiltMax] <- 0

  NormExp <- NormExp[rowSums(NormExp==0) < ncol(NormExp), ]

  pheatmap::pheatmap(NormExp.Bin*1, cluster_cols = FALSE)

  NormExp.Bin <- NormExp > apply(NormExp, 1, min) + (apply(NormExp, 1, max) - apply(NormExp, 1, min))*Thr

  pheatmap::pheatmap(1*NormExp.Bin[order(rowMeans(t(t(NormExp.Bin)*1:ncol(NormExp.Bin)))), ],
                     cluster_cols = FALSE, cluster_rows = FALSE)

  SinglePk <- lapply(apply(NormExp.Bin, 1, which), function(x){
    if(all(min(x):max(x) %in% x)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })

  PlotMat <- 1*NormExp.Bin[unlist(SinglePk),]

  PlotMat.Numb <- t(t(PlotMat)*1:ncol(PlotMat))
  PlotMat.Numb[PlotMat.Numb == 0] <- NA

  pheatmap::pheatmap(PlotMat[order(apply(PlotMat.Numb, 1, median, na.rm=TRUE)),], show_colnames = FALSE,
                     show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)

  barplot(colSums(NormExp.Bin[unlist(SinglePk),]))

  StageGenes <- lapply(ProcStruct$StageOnNodes, function(x){
    # SelX <- intersect(paste(x), colnames(GeneExprMat.DF.Split.Mean.Mat.Bin))
    SelX <- intersect(x, 1:ncol(NormExp.Bin))
    if(length(SelX)>0){
      if(length(SelX)==1){
        NormExp.Bin[unlist(SinglePk), SelX]
      } else {
        apply(NormExp.Bin[unlist(SinglePk), SelX], 1, function(y){
          any(which(y))
        })
      }
    } else {
      NA
    }

  })

  StageGenes.Names <- lapply(StageGenes, function(x){
    names(which(x))
  })

  barplot(unlist(lapply(StageGenes.Names, length)), las = 2, horiz = FALSE)

  return(StageGenes.Names)

}


















#
#
#
# GetGenesWithPeaks <- function(DataStruct, ProcStruct, FiltMax = 1, Thr = 3) {
#
#   # WorkPG <- DataStruct$Analysis$PGStructs[[length(DataStruct$Analysis$PGStructs)]]
#   # TaxList <- WorkPG$TaxonList[[length(WorkPG$TaxonList)]]
#   # TaxVect <- rep(NA, nrow(WorkPG$Data))
#   # names(TaxVect) <- rownames(WorkPG$Data)
#   #
#   # lapply(1:length(TaxList), function(j){
#   #   TaxVect[ TaxList[[j]] ] <<- j
#   # })
#   #
#   #
#   # GeneExprMat.DF <- data.frame(t(DataStruct$ExpMat[, names(ProcStruct$CellsPT)]))
#   # colnames(GeneExprMat.DF) <- rownames(DataStruct$ExpMat)
#   #
#   # GeneExprMat.DF.Split <- split(GeneExprMat.DF[names(TaxVect), ], TaxVect[names(TaxVect)])
#   #
#   # GeneExprMat.DF.Split.Mean <- lapply(GeneExprMat.DF.Split, colMeans)
#   # GeneExprMat.DF.Split.Mean.Mat <- NULL
#   #
#   # for(i in 1:length(GeneExprMat.DF.Split.Mean)){
#   #   GeneExprMat.DF.Split.Mean.Mat <- cbind(GeneExprMat.DF.Split.Mean.Mat, GeneExprMat.DF.Split.Mean[[i]])
#   # }
#   #
#   # GeneExprMat.DF.Split.Mean.Mat <- GeneExprMat.DF.Split.Mean.Mat[rowSums(GeneExprMat.DF.Split.Mean.Mat) > 0, ]
#   #
#   # colnames(GeneExprMat.DF.Split.Mean.Mat) <- paste((1:length(TaxList))[unlist(lapply(TaxList, function(x){!any(is.na(x))}))])
#   #
#   # RefPath <- ProcStruct$ExtPath[-length(ProcStruct$ExtPath)]
#   # ReordCol <- RefPath[RefPath %in% colnames(GeneExprMat.DF.Split.Mean.Mat)]
#   #
#   # GeneExprMat.DF.Split.Mean.Mat <- GeneExprMat.DF.Split.Mean.Mat[, ReordCol]
#   #
#   # GeneMin <- apply(GeneExprMat.DF.Split.Mean.Mat, 1, min)
#   # GeneMax <- apply(GeneExprMat.DF.Split.Mean.Mat, 1, max)
#   # GeneMean <- apply(GeneExprMat.DF.Split.Mean.Mat, 1, mean)
#   # GeneSd <- apply(GeneExprMat.DF.Split.Mean.Mat, 1, sd)
#   # GeneMedian <- apply(GeneExprMat.DF.Split.Mean.Mat, 1, median)
#   #
#   # GeneExprMat.DF.Split.Mean.Mat <- GeneExprMat.DF.Split.Mean.Mat[(GeneMax - GeneMin) > FiltMax, ]
#
#   # ToRem <- which(GeneMAD == 0)
#
#   # MadMed <- (GeneExprMat.DF.Split.Mean.Mat[-ToRem, ] - GeneMedian[-ToRem])/GeneMAD[-ToRem]
#
#   # GeneExprMat.DF.Split.Mean.Mat <- GeneExprMat.DF.Split.Mean.Mat - apply(GeneExprMat.DF.Split.Mean.Mat, 1, min)
#   # GeneExprMat.DF.Split.Mean.Mat <- GeneExprMat.DF.Split.Mean.Mat/apply(GeneExprMat.DF.Split.Mean.Mat, 1, max)
#
#   NormExp <- ProcStruct$NodesExp
#   NormExp[NormExp < .5] <- 0
#
#   NormExp <- NormExp[rowSums(NormExp==0) < ncol(NormExp), ]
#
#   pheatmap::pheatmap(t(NormExp.Bin*1), cluster_cols = FALSE)
#
#   NormExp.Bin <- NormExp > apply(NormExp, 1, quantile, .8)
#   pheatmap::pheatmap(1*NormExp.Bin[order(rowMeans(t(t(NormExp.Bin)*1:ncol(NormExp.Bin)))), ],
#                      cluster_cols = FALSE, cluster_rows = FALSE)
#
#
#
#
#   -NormExp <- (ProcStruct$NodesExp - apply(ProcStruct$NodesExp, 1, min))
#
#   NormExp <- NormExp/apply(NormExp, 1, max)
#
#   GeneExprMat.DF.Split.Mean.Mat.Bin <- NormExp > Thr
#
#   # GeneExprMat.DF.Split.Mean.Mat.Bin <- GeneExprMat.DF.Split.Mean.Mat > Thr
#
#   # GeneExprMat.DF.Split.Mean.Mat.Bin <- MadMed > Thr
#
#   GeneExprMat.DF.Split.Mean.Mat.Bin <- GeneExprMat.DF.Split.Mean.Mat.Bin[apply(GeneExprMat.DF.Split.Mean.Mat.Bin, 1, any), ]
#
#   SinglePk <- lapply(apply(GeneExprMat.DF.Split.Mean.Mat.Bin, 1, which), function(x){
#     if(all(min(x):max(x) %in% x)){
#       return(TRUE)
#     } else {
#       return(FALSE)
#     }
#   })
#
#   PlotMat <- 1*GeneExprMat.DF.Split.Mean.Mat.Bin[unlist(SinglePk),]
#
#   PlotMat.Numb <- t(t(PlotMat)*1:ncol(PlotMat))
#   PlotMat.Numb[PlotMat.Numb == 0] <- NA
#
#   pheatmap::pheatmap(PlotMat[order(apply(PlotMat.Numb, 1, median, na.rm=TRUE)),], show_colnames = FALSE,
#                      show_rownames = FALSE, cluster_cols = FALSE, cluster_rows = FALSE)
#
#   barplot(colSums(GeneExprMat.DF.Split.Mean.Mat.Bin[unlist(SinglePk),]))
#
#   StageGenes <- lapply(ProcStruct$StageOnNodes, function(x){
#     # SelX <- intersect(paste(x), colnames(GeneExprMat.DF.Split.Mean.Mat.Bin))
#     SelX <- intersect(x, 1:ncol(GeneExprMat.DF.Split.Mean.Mat.Bin))
#     if(length(SelX)>0){
#       if(length(SelX)==1){
#         GeneExprMat.DF.Split.Mean.Mat.Bin[unlist(SinglePk), SelX]
#       } else {
#         apply(GeneExprMat.DF.Split.Mean.Mat.Bin[unlist(SinglePk), SelX], 1, function(y){
#           any(which(y))
#         })
#       }
#     } else {
#       NA
#     }
#
#   })
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


















PeackGenes.Kowa <-  GetGenesWithPeaks(DataStruct = Data.Kowa,
                                      ProcStruct = Proc.Exp.Kowa,
                                      FiltMax = .5, Thr = .9)
PeackGenes.Buet <-  GetGenesWithPeaks(DataStruct = Data.Buet,
                                      ProcStruct = Proc.Exp.Buet,
                                      FiltMax = .5, Thr = .9)
PeackGenes.Sasa <-  GetGenesWithPeaks(DataStruct = Data.Sasa,
                                      ProcStruct = Proc.Exp.Sasa,
                                      FiltMax = .5, Thr = .9)














names(PeackGenes.Kowa)
names(PeackGenes.Buet)
names(PeackGenes.Sasa)


G0.kowa <- unlist(PeackGenes.Kowa[c(1,10)])


G1.kowa <- unlist(PeackGenes.Kowa[2:6])
G1.buet <- unlist(PeackGenes.Buet[1:2])
G1.sasa <- unlist(PeackGenes.Sasa[1:2])

S.kowa <- unlist(PeackGenes.Kowa[6:8])
S.buet <- unlist(PeackGenes.Buet[2:4])
S.sasa <- unlist(PeackGenes.Sasa[2:4])

G2M.kowa <- unlist(PeackGenes.Kowa[8:10])
G2M.buet <- unlist(PeackGenes.Buet[4:6])
G2M.sasa <- unlist(PeackGenes.Sasa[4:6])


G0.All <- G0.kowa

G1.All <- unique(c(G1.kowa, G1.buet, G1.sasa))

S.All <- c(S.kowa, S.buet, S.sasa)

G2.All <- c(G2M.kowa, G2M.buet, G2M.sasa)


ToRemove <- c(intersect(G1.All, S.All),
              intersect(S.All, G2.All),
              intersect(G1.All, G2.All))


G1.Sel <- G1.All
G1.Sel <- G1.Sel[!(G1.Sel %in% ToRemove)]

S.Sel <- S.All
S.Sel <- S.Sel[!(S.Sel %in% ToRemove)]

G2.Sel <- G2.All
G2.Sel <- G2.Sel[!(G2.Sel %in% ToRemove)]

G0.Sel <- G0.All
G0.Sel <- G0.Sel[!(G0.Sel %in% ToRemove)]


length(G1.Sel)
length(S.Sel)
length(G2.Sel)
length(G0.Sel)


AllMat.All <- cbind(AllMat.Sasa, AllMat.Buet)

AllMat.All <- AllMat.All[, apply(AllMat.All, 2, any)]

pheatmap::pheatmap(1*unique(AllMat.All))

pheatmap::pheatmap(1*(AllMat.All))
























# pheatmap::pheatmap(cor(Proc.Exp.Buet$NodesExp, method = "spe"))

# pheatmap::pheatmap(cor(t(Proc.Exp.Buet$NodesExp), method = "spe"), show_rownames = FALSE, show_colnames = FALSE)


Max.Bue <- apply(Proc.Exp.Buet$NodesExp, 1, which.max)
Max.Sasa <- apply(Proc.Exp.Sasa$NodesExp, 1, which.max)

Shared <- intersect(names(Max.Bue), names(Max.Sasa))

plot(Max.Bue[Shared], Max.Sasa[Shared])

cor.test(Max.Bue[Shared], Max.Sasa[Shared], method = "pea")

AllGenes <- unique(c(unlist(PeackGenes.Kowa), unlist(PeackGenes.Buet), unlist(PeackGenes.Sasa)))

AllMat.Kowa <- sapply(PeackGenes.Kowa, function(x){AllGenes %in% x})
colnames(AllMat.Kowa) <- paste(colnames(AllMat.Kowa), "_Kowa", sep='')

AllMat.Buet <- sapply(PeackGenes.Buet, function(x){AllGenes %in% x})
colnames(AllMat.Buet) <- paste(colnames(AllMat.Buet), "_Buet", sep='')

AllMat.Sasa <- sapply(PeackGenes.Sasa, function(x){AllGenes %in% x})
colnames(AllMat.Sasa) <- paste(colnames(AllMat.Sasa), "_Sasa", sep='')

# AllMat.All <- cbind(AllMat.Kowa, AllMat.Buet, AllMat.Sasa)

AllMat.All <- cbind(AllMat.Sasa, AllMat.Buet)

AllMat.All <- AllMat.All[, apply(AllMat.All, 2, any)]

pheatmap::pheatmap(1*unique(AllMat.All))

pheatmap::pheatmap(1*(AllMat.All))


IntMat <- apply(AllMat.All, 2, function(x){
  apply(AllMat.All, 2, function(y){ sum(x & y)/min(sum(x), sum(y)) })
})

diag(IntMat) <- NA

pheatmap::pheatmap(IntMat[rownames(IntMat) %in% c("G1_Buet", "G1+S_Buet"),
                          !(colnames(IntMat) %in% c("G1_Buet", "G1+S_Buet"))],
                   cluster_rows = TRUE, cluster_cols = TRUE)


pheatmap::pheatmap(IntMat[rownames(IntMat) %in% c("G1+S_Buet", "S_Buet", "S+G2M_Buet"),
                          !(colnames(IntMat) %in% c("G1+S_Buet", "S_Buet", "S+G2M_Buet"))],
                   cluster_rows = TRUE, cluster_cols = TRUE)


pheatmap::pheatmap(IntMat[rownames(IntMat) %in% c("G2M_Buet", "S+G2M_Buet"),
                          !(colnames(IntMat) %in% c("G2M_Buet", "S+G2M_Buet"))],
                   cluster_rows = TRUE, cluster_cols = TRUE)








AllMat.All <- cbind(AllMat.Kowa, AllMat.Buet)

AllMat.All <- AllMat.All[, apply(AllMat.All, 2, any)]

IntMat <- apply(AllMat.All, 2, function(x){
  apply(AllMat.All, 2, function(y){ sum(x & y)/min(sum(x), sum(y)) })
})

diag(IntMat) <- NA

pheatmap::pheatmap(IntMat[rownames(IntMat) %in% c("G1_Buet", "G1+S_Buet"),
                          !(colnames(IntMat) %in% c("G1_Buet", "G1+S_Buet"))],
                   cluster_rows = TRUE, cluster_cols = TRUE)


pheatmap::pheatmap(IntMat[rownames(IntMat) %in% c("G1+S_Buet", "S_Buet", "S+G2M_Buet"),
                          !(colnames(IntMat) %in% c("G1+S_Buet", "S_Buet", "S+G2M_Buet"))],
                   cluster_rows = TRUE, cluster_cols = TRUE)


pheatmap::pheatmap(IntMat[rownames(IntMat) %in% c("G2M_Buet", "S+G2M_Buet"),
                          !(colnames(IntMat) %in% c("G2M_Buet", "S+G2M_Buet"))],
                   cluster_rows = TRUE, cluster_cols = TRUE)





















diag(IntMat) <- 0

rownames(IntMat)[apply(IntMat, 1, which.max)]



pheatmap::pheatmap(1*(IntMat > .3),
                   clustering_distance_rows = "binary",
                   clustering_distance_cols = "binary")













FilterRepeating <- function(PKG.List) {

  PKG.All <- unique(unlist(PKG.List))

  GeneMat.List <- lapply(PKG.List, function(x) {
    PKG.All %in% x
  })

  GeneMat <- NULL

  for(i in 1:length(GeneMat.List)){
    GeneMat <- rbind(GeneMat, GeneMat.List[[i]])
  }

  colnames(GeneMat) <- PKG.All

  MultPh <- which(colSums(GeneMat)>1)

  if(length(MultPh) == 0){
    return(PKG.List)
  }

  print(MultPh)

  MultStages <- lapply(apply(GeneMat[,MultPh], 2, which), function(x) {
    names(PKG.List)[x]
  })

  DividedStages <- lapply(MultStages, function(x) {
    strsplit(x, "+", fixed = TRUE)
  })

  SingleStageGenes <- lapply(DividedStages, function(x){
    StagesAll <- rep(NA, length(unique(unlist(x))))
    names(StagesAll) <- unique(unlist(x))

    for(Stg in names(StagesAll)){
      ListBel <- lapply(x, function(x){
        Stg %in% x
      })

      StagesAll[Stg] <- all(unlist(ListBel))

    }

    if(any(StagesAll)){
      return(which(StagesAll))
    } else {
      return(NA)
    }

  })

  PKG.List <- lapply(PKG.List, function(x){
    x[!(x %in% names(MultPh))]
  })

  SingleStageGenes <- SingleStageGenes[!unlist(lapply(SingleStageGenes, function(x){any(is.na(x))}))]

  if(length(SingleStageGenes) > 0){
    for(i in 1:length(SingleStageGenes)){
      PKG.List[[ SingleStageGenes[[i]] ]] <- union(PKG.List[[ SingleStageGenes[[i]] ]],
                                                   names(SingleStageGenes)[i])
    }
  }

  return(PKG.List)
}







PeackGenes.Kowa <-  FilterRepeating(PeackGenes.Kowa)
PeackGenes.Buet <-  FilterRepeating(PeackGenes.Buet)
PeackGenes.Sasa <-  FilterRepeating(PeackGenes.Sasa)


barplot(unlist(lapply(PeackGenes.Kowa, length)), las = 2)
barplot(unlist(lapply(PeackGenes.Buet, length)), las = 2)
barplot(unlist(lapply(PeackGenes.Sasa, length)), las = 2)



AllG1 <- c(PeackGenes.Kowa$`G1(late)`,
           PeackGenes.Kowa$`G1(early)`,
           PeackGenes.Kowa$`G1(early)+G1(late)`,
           PeackGenes.Kowa$`G1(late)`,
           PeackGenes.Kowa$`G1(late)+S`,
           PeackGenes.Buet$G1,
           PeackGenes.Buet$`G1+S`,
           PeackGenes.Sasa$G1,
           PeackGenes.Sasa$`G1+S`)

G1Mat <- rbind(AllG1 %in% c(PeackGenes.Kowa$`G1(late)`,
                            PeackGenes.Kowa$`G1(early)`,
                            PeackGenes.Kowa$`G1(early)+G1(late)`,
                            PeackGenes.Kowa$`G1(late)`,
                            PeackGenes.Kowa$`G1(late)+S`),
               AllG1 %in% c(PeackGenes.Buet$G1,
                            PeackGenes.Buet$`G1+S`),
               AllG1 %in% c(PeackGenes.Sasa$G1,
                            PeackGenes.Sasa$`G1+S`))

colnames(G1Mat) <- AllG1

table(colSums(G1Mat))




AllS <- c(PeackGenes.Kowa$S,
          PeackGenes.Buet$S,
          PeackGenes.Sasa$S)

SMat <- rbind(AllS %in% PeackGenes.Kowa$S,
              AllS %in% PeackGenes.Buet$S,
              AllS %in% PeackGenes.Sasa$S)

colnames(SMat) <- AllS

table(colSums(SMat))


AllG2M <- c(PeackGenes.Kowa$`G2/M`,
          PeackGenes.Buet$G2M,
          PeackGenes.Sasa$G2M)

G2MMat <- rbind(AllG2M %in% PeackGenes.Kowa$`G2/M`,
                AllG2M %in% PeackGenes.Buet$G2M,
                AllG2M %in% PeackGenes.Sasa$G2M)

colnames(G2MMat) <- AllG2M

table(colSums(G2MMat))


AllGenes <- unique(c(unlist(PeackGenes.Kowa), unlist(PeackGenes.Buet), unlist(PeackGenes.Sasa)))

AllMat.Kowa <- sapply(PeackGenes.Kowa, function(x){AllGenes %in% x})
colnames(AllMat.Kowa) <- paste(colnames(AllMat.Kowa), "_Kowa", sep='')

AllMat.Buet <- sapply(PeackGenes.Buet, function(x){AllGenes %in% x})
colnames(AllMat.Buet) <- paste(colnames(AllMat.Buet), "_Buet", sep='')

AllMat.Sasa <- sapply(PeackGenes.Sasa, function(x){AllGenes %in% x})
colnames(AllMat.Sasa) <- paste(colnames(AllMat.Sasa), "_Sasa", sep='')

AllMat.All <- cbind(AllMat.Kowa, AllMat.Buet, AllMat.Sasa)

AllMat.All <- AllMat.All[, apply(AllMat.All, 2, any)]

pheatmap::pheatmap(1*unique(AllMat.All))

pheatmap::pheatmap(1*(AllMat.All))


IntMat <- apply(AllMat.All, 2, function(x){
  apply(AllMat.All, 2, function(y){ sum(x & y)/min(sum(x), sum(y)) })
})

diag(IntMat) <- NA

pheatmap::pheatmap(1*(IntMat>.5))


# AllFreman_Mouse
#
#
# Genes.Buett <- Data.Buet$Analysis$Genes[[length(Data.Buet$Analysis$Genes)]]
# Genes.Sasa <- Data.Sasa$Analysis$Genes[[length(Data.Sasa$Analysis$Genes)]]
# Genes.Kowa <- Data.Kowa$Analysis$Genes[[length(Data.Kowa$Analysis$Genes)]]
#
# length(Genes.Buett)
# length(Genes.Kowa)
# length(Genes.Sasa)
#
# library(VennDiagram)
#
# grid.newpage()
# draw.quad.venn(area1 = length(Genes.Buett), area2 = length(Genes.Kowa),
#                area3 = length(Genes.Sasa), area4 = length(MouseGenes_GOCellCycle),
#                n12 = length(intersect(Genes.Buett, Genes.Kowa)),
#                n23 = length(intersect(Genes.Sasa, Genes.Kowa)),
#                n13 = length(intersect(Genes.Buett, Genes.Sasa)),
#                n14 = length(intersect(Genes.Buett, MouseGenes_GOCellCycle)),
#                n24 = length(intersect(MouseGenes_GOCellCycle, Genes.Kowa)),
#                n34 = length(intersect(Genes.Sasa, MouseGenes_GOCellCycle)),
#                n123 = length(intersect(Genes.Buett, intersect(Genes.Kowa, Genes.Sasa))),
#                n124 = length(intersect(Genes.Buett, intersect(Genes.Kowa, MouseGenes_GOCellCycle))),
#                n134 = length(intersect(Genes.Buett, intersect(Genes.Sasa, MouseGenes_GOCellCycle))),
#                n234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, MouseGenes_GOCellCycle))),
#                n1234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, intersect(Genes.Buett, MouseGenes_GOCellCycle)))),
#                category = c("Buettner et al.", "Kowalczyk et al.", "Sasagawa et al.", "GO"),
#                fill = c("red", "green", "blue", "white"), cex = 3)
# grid.newpage()
#
#
#
#
#
# Freeman_G1S_CC4
#
#
#
#
#
#
