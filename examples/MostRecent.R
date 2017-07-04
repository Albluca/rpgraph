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
  AllTerms <- sort(c("GO:0007049", AllTerms))

  MouseGenes_GOCellCycle <- getBM(attributes = c("external_gene_name"),
                                  filters = "go_id",
                                  values = AllTerms,
                                  mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))

  MouseGenes_GOCellCycle <- unlist(MouseGenes_GOCellCycle, use.names = FALSE)

  AllWit_Mouse <- gsub("(\\b[a-z]{1})", "\\U\\1", tolower(unlist(StageAssociation_Whit_Ext[-c(1,2)], use.names = FALSE)), perl=TRUE)

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

  AllWit_Mouse <- gsub("(\\b[a-z]{1})", "\\U\\1", tolower(unlist(StageAssociation_Whit_Ext[-c(1,2)], use.names = FALSE)), perl=TRUE)

}

# Define functions ----------------------------------------------------------------














# Kowalczyk et al ----------------------------------------------------------------

# Load data

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




# SelectedSamples <- Murine_HSC_Samples$`cell ID`[Murine_HSC_Samples$population %in% c("'old ST-HSC'", "'old ST-HSC(rep)'", "'young ST-HSC'")]



SelectedSamples <- intersect(SelectedSamples, colnames(Murine_HSC_Data))

Categories <- CCVect[paste(SelectedSamples)]

AllData <- data.matrix(Murine_HSC_Data[-1,paste(SelectedSamples)])
rownames(AllData) <- gsub("'", '', unlist(Murine_HSC_Data[-1, 1]))

GeneExprMat <- data.matrix(AllData)
StartSet <- MouseGenes_GOCellCycle
Categories <- Categories


FullExpData.Kowa <- log10(GeneExprMat + 1)
FullCat.Kowa <- Categories






# GeneExprMat, StartSet, Categories, Topology = "Circle", IgnoreTail = FALSE, PlanVarLimit = .9, PlanVarLimitIC = .95,
# DistillThr = .7, Log = TRUE, StartQuant = .5, OutThr = 3, OutThrPCA = 3, Title = '', MinProlCells = 20, PCACenter = FALSE,
# PlotDebug = FALSE, Mode = "VarPC", ExtMode = 2, nNodes = 40, InitStructNodes = 20, DipPVThr = 1e-4, PCAFilter = TRUE,
# PCAProjCenter = TRUE, StopCrit = .95, Filter = TRUE, QuaThr = .5, CompareSet = list(), EstProlif = "MeanPerc",
# KeepOrgProlifiration = FALSE, NonG0Cell = NULL







Output.Kowa <- SelectGenesOnGraph(DataSet = GeneExprMat, StartSet = MouseGenes_GOCellCycle, Categories = Categories,
                             PCACenter = TRUE, PCAProjCenter = TRUE)


write_rds(x = list(Analysis = Output.Kowa,
                   ExpMat = log10(GeneExprMat + 1),
                   Cats = Categories),
          path = "~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.rds")






TargetStruct <- Output$PGStructs[[length(Output$PGStructs)]]

Proc.Exp <- PlotOnStages(Structure = "Circle",
                         Categories = TargetStruct$Categories,
                         nGenes = 10,
                         TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                         PrinGraph = TargetStruct$PrinGraph,
                         Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                         SelThr = .32,
                         ComputeOverlaps = TRUE,
                         ExpData = TargetStruct$FiltExp,
                         RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                         PCACenter = TargetStruct$PCAData$center,
                         PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                         OrderOnCat = TRUE,
                         SmoothPoints = 2, MinCellPerNode = 6)


NormGE <- apply(Proc.Exp$NodesExp, 1, function(x){
  x1 <- x-min(x)
  x1 <- x1/max(x1)
  return(x1)
})

NormGE <- NormGE > .95

Peacks <- apply(NormGE, 2, which)

PeackStages <- lapply(Peacks, function(x){
  sapply(Proc.Exp$StageOnNodes, function(y){
    any(x %in% y)
  })
})

MatExp <- matrix(unlist(PeackStages), nrow = length(Proc.Exp$StageOnNodes))
colnames(MatExp) <- names(PeackStages)
rownames(MatExp) <- names(Proc.Exp$StageOnNodes)

Peacks <- apply(MatExp, 1, which)




pheatmap::pheatmap(t(NormGE), cluster_cols = FALSE)

for(i in 1:length(Proc.Exp$StageOnNodes)){
  if(length(unlist(Proc.Exp$StageOnNodes[i]))>0){
    pheatmap::pheatmap(t(NormGE[unlist(Proc.Exp$StageOnNodes[i]), ]),
                       cluster_cols = FALSE, main = names(Proc.Exp$StageOnNodes)[i])
  }
}




#
#
# Structure = "Circle"
# Categories = TargetStruct$Categories
# nGenes = 2
# TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]]
# PrinGraph = Output$PGStructs[[length(Output$PGStructs)]]$PrinGraph
# Net = TargetStruct$Net[[length(TargetStruct$Net)]]
# SelThr = .39
# ComputeOverlaps = FALSE
# ExpData = TargetStruct$FiltExp
# RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims]
# PCACenter = TargetStruct$PCAData$center
# PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]]
# OrderOnCat = TRUE
# SmoothPoints = 2
#







# SampledList <- lapply(as.list(1:nRep), function(i){
#   SampVect <- sample(1:length(TaxList), size = length(TaxVect), replace = TRUE)
#   names(SampVect) <- names(TaxVect)
#   return(SampVect)
# })
#
#
# cl = makeCluster(4)
# clusterExport(cl, c("GeneExprMat.DF", "Steps"))
# RandData <- parLapply(cl, SampledList, DoStuff)
# stopCluster(cl)
#
#
# RandDataMat <- NULL
#
# for(i in 1:length(RandData)){
#   RandDataMat <- rbind(RandDataMat, RandData[[i]])
# }
#
#
# AllGPV <- unlist(lapply(apply(t(RandDataMat)-Base, 1, wilcox.test, alternative = "greater"), "[[", "p.value"))
# hist(AllGPV)

#
# AddGeneNames <- names(sort(AllGPV[AllGPV > .05], decreasing = TRUE))[1:round(length(UsedGenes[[2]])/10)]
#
# length(setdiff(AddGeneNames, UsedGenes[[2]]))
# length(AddGeneNames)
#
#
#
# UsedGenes[[3]] <- union(AddGeneNames, UsedGenes[[2]])
#
#
# Steps[[3]] <- ProjectAndCompute(DataSet = GeneExprMat, GeneSet = UsedGenes[[3]], VarThr = .99, nNodes = 40, Log = TRUE, Categories = Categories,
#                                               Filter = TRUE, OutThr = 5, PCAFilter = TRUE, OutThrPCA = 3, GraphType = "Circle", PlanVarLimit = .85, PlanVarLimitIC = .9,
#                                               MinBranDiff = 3, InitStructNodes = 20, ForceLasso = FALSE, EstProlif = "MeanPerc", QuaThr = .5, NonG0Cell = Steps[[1]]$NonG0Cell, DipPVThr = 1e-4,
#                                               MinProlCells = 50, PCACenter = TRUE, PCAProjCenter = TRUE, PlotDebug = FALSE, PlotIntermediate = FALSE)
#
#
#
