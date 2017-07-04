library(readr)
library(Pricecycle)
library(GO.db)
library(biomaRt)

AllTerms <- GOBPOFFSPRING[["GO:0007049"]]
AllTerms <- sort(c("GO:0007049", AllTerms))

MouseGenes_GOCellCycle <- getBM(attributes = c("external_gene_name"),
                                filters = "go",
                                values = AllTerms,
                                mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"))

MouseGenes_GOCellCycle <- MouseGenes_GOCellCycle[,1]





Human_GOCellCycle <- getBM(attributes = c("external_gene_name"),
                           filters = "go",
                           values = AllTerms,
                           mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"))

Human_GOCellCycle <- Human_GOCellCycle[,1]






# Load data

Data.Sasa <- read_rds("~/Google Drive/Datasets/Sasagawa et al - Murine Stem Cells (EB5 cell line)/rPG_Last.rds")
Data.Kowa <- read_rds("~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last.rds")
Data.Buet <- read_rds("~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/rPG_Last.rds")

Data.Sasa.Norm <- read_rds("~/Google Drive/Datasets/Sasagawa et al - Murine Stem Cells (EB5 cell line)/rPG_Last-Filtered.rds")
Data.Kowa.Norm <- read_rds("~/Google Drive/Datasets/Kowalczyk et al - Murine hematopoietic stem cells/rPG_Last-Filtered.rds")
Data.Buet.Norm <- read_rds("~/Google Drive/Datasets/Buettner et al. - Murie embrionic stem cells/rPG_Last-Filtered.rds")

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
                              SmoothPoints = 2,
                              MinCellPerNode = 2,
                              Title = 'Kowalczyk et al')

ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Kowalczyk et al.",
  ExpValues = NULL, PlotPC3 = FALSE,
  PointAlpha = .5
)


ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Kowalczyk et al.",
  ExpValues = NULL, PlotPC3 = TRUE,
  PointAlpha = .5
)




TargetStruct <- Data.Buet$Analysis$PGStructs[[length(Data.Buet$Analysis$PGStructs)]]
Proc.Exp.Buet <- PlotOnStages(Structure = "Circle",
                              Categories = TargetStruct$Categories,
                              nGenes = 2,
                              TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                              PrinGraph = TargetStruct$PrinGraph,
                              Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                              SelThr = .34,
                              ComputeOverlaps = TRUE,
                              ExpData = TargetStruct$FiltExp,
                              RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                              PCACenter = TargetStruct$PCAData$center,
                              PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                              OrderOnCat = TRUE,
                              SmoothPoints = 2,
                              MinCellPerNode = 2,
                              Title = 'Buettner et al')




ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Buettner et al.",
  ExpValues = NULL, PlotPC3 = FALSE,
  PointAlpha = .5
)


ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Buettner et al.",
  ExpValues = NULL, PlotPC3 = TRUE,
  PointAlpha = .5
)






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
                              SmoothPoints = 0,
                              MinCellPerNode = 1,
                              Title = 'Sasagawa et al')




ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Sasagawa et al.",
  ExpValues = NULL, PlotPC3 = FALSE,
  PointAlpha = .5
)


ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Sasagawa et al.",
  ExpValues = NULL, PlotPC3 = TRUE,
  PointAlpha = .5
)














TargetStruct <- Data.Kowa.Norm$Analysis$PGStructs[[length(Data.Kowa.Norm$Analysis$PGStructs)]]
Proc.Exp.Kowa.Norm <- PlotOnStages(Structure = "Circle",
                              Categories = TargetStruct$Categories,
                              nGenes = 2,
                              TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                              PrinGraph = TargetStruct$PrinGraph,
                              Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                              SelThr = .292,
                              ComputeOverlaps = TRUE,
                              ExpData = TargetStruct$FiltExp,
                              RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                              PCACenter = TargetStruct$PCAData$center,
                              PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                              OrderOnCat = TRUE,
                              SmoothPoints = 1,
                              MinCellPerNode = 2,
                              Title = 'Kowalczyk et al')

ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Kowalczyk et al. (Norm)",
  ExpValues = NULL, PlotPC3 = FALSE,
  PointAlpha = .5
)


ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Kowalczyk et al. (Norm)",
  ExpValues = NULL, PlotPC3 = TRUE,
  PointAlpha = .5
)




TargetStruct <- Data.Buet.Norm$Analysis$PGStructs[[length(Data.Buet.Norm$Analysis$PGStructs)]]
Proc.Exp.Buet.Norm <- PlotOnStages(Structure = "Circle",
                              Categories = TargetStruct$Categories,
                              nGenes = 2,
                              TaxonList = TargetStruct$TaxonList[[length(TargetStruct$TaxonList)]],
                              PrinGraph = TargetStruct$PrinGraph,
                              Net = TargetStruct$Net[[length(TargetStruct$Net)]],
                              SelThr = .34,
                              ComputeOverlaps = TRUE,
                              ExpData = TargetStruct$FiltExp,
                              RotatioMatrix = TargetStruct$PCAData$rotation[,1:TargetStruct$nDims],
                              PCACenter = TargetStruct$PCAData$center,
                              PointProjections = TargetStruct$ProjPoints[[length(TargetStruct$ProjPoints)]],
                              OrderOnCat = TRUE,
                              SmoothPoints = 2,
                              MinCellPerNode = 2,
                              Title = 'Buettner et al')




ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Buettner et al. (Norm)",
  ExpValues = NULL, PlotPC3 = FALSE,
  PointAlpha = .5
)


ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Buettner et al. (Norm)",
  ExpValues = NULL, PlotPC3 = TRUE,
  PointAlpha = .5
)






TargetStruct <- Data.Sasa.Norm$Analysis$PGStructs[[length(Data.Sasa.Norm$Analysis$PGStructs)]]
Proc.Exp.Sasa.Norm <- PlotOnStages(Structure = "Circle",
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
                              SmoothPoints = 0,
                              MinCellPerNode = 1,
                              Title = 'Sasagawa et al')




ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Sasagawa et al. (Norm)",
  ExpValues = NULL, PlotPC3 = FALSE,
  PointAlpha = .5
)


ProjectOnCircle(
  Points = TargetStruct$Data,
  Edges = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Edges,
  Nodes = TargetStruct$IntGrahs[[length(TargetStruct$IntGrahs)]]$Nodes,
  Categories = TargetStruct$Categories,
  Title = "Sasagawa et al. (Norm)",
  ExpValues = NULL, PlotPC3 = TRUE,
  PointAlpha = .5
)
























# Create structure to use thge dashboard


InputList <- list(
  list(Name = "Buettner et al", Expression = Data.Buet$ExpMat,
                       Categories = Data.Buet$Cats,
                       OrderedData = Proc.Exp.Buet,
                       PGStruct = Data.Buet$Analysis$PGStructs[[length(Data.Buet$Analysis$PGStructs)]],
                       TaxonList = Data.Buet$Analysis$PGStructs[[length(Data.Buet$Analysis$PGStructs)]]$TaxonList
       ),
  list(Name = "Buettner et al (Norm)", Expression = Data.Buet.Norm$ExpMat,
       Categories = Data.Buet.Norm$Cats,
       OrderedData = Proc.Exp.Buet.Norm,
       PGStruct = Data.Buet.Norm$Analysis$PGStructs[[length(Data.Buet.Norm$Analysis$PGStructs)]],
       TaxonList = Data.Buet.Norm$Analysis$PGStructs[[length(Data.Buet.Norm$Analysis$PGStructs)]]$TaxonList
  ),
  list(Name = "Kowalczyk et al", Expression = Data.Kowa$ExpMat,
       Categories = Data.Kowa$Cats,
       OrderedData = Proc.Exp.Kowa,
       PGStruct = Data.Kowa$Analysis$PGStructs[[length(Data.Kowa$Analysis$PGStructs)]],
       TaxonList = Data.Kowa$Analysis$PGStructs[[length(Data.Kowa$Analysis$PGStructs)]]$TaxonList
  ),
list(Name = "Kowalczyk et al (Norm)", Expression = Data.Kowa.Norm$ExpMat,
     Categories = Data.Kowa.Norm$Cats,
     OrderedData = Proc.Exp.Kowa.Norm,
     PGStruct = Data.Kowa.Norm$Analysis$PGStructs[[length(Data.Kowa.Norm$Analysis$PGStructs)]],
     TaxonList = Data.Kowa.Norm$Analysis$PGStructs[[length(Data.Kowa.Norm$Analysis$PGStructs)]]$TaxonList
),
list(Name = "Sasagawa et al", Expression = Data.Sasa$ExpMat,
     Categories = Data.Sasa$Cats,
     OrderedData = Proc.Exp.Sasa,
     PGStruct = Data.Sasa$Analysis$PGStructs[[length(Data.Sasa$Analysis$PGStructs)]],
     TaxonList = Data.Sasa$Analysis$PGStructs[[length(Data.Sasa$Analysis$PGStructs)]]$TaxonList
),
list(Name = "Sasagawa et al (Norm)", Expression = Data.Sasa.Norm$ExpMat,
     Categories = Data.Sasa.Norm$Cats,
     OrderedData = Proc.Exp.Sasa.Norm,
     PGStruct = Data.Sasa.Norm$Analysis$PGStructs[[length(Data.Sasa.Norm$Analysis$PGStructs)]],
     TaxonList = Data.Sasa.Norm$Analysis$PGStructs[[length(Data.Sasa.Norm$Analysis$PGStructs)]]$TaxonList
)
)


SpanVect <- c(.1, .1, .05, .05, .3, .3)


for(i in 1:length(InputList)){

  p <- PlotOnPseudotime(WorkStruct = InputList[[i]]$OrderedData,
                        Expression = InputList[[i]]$Expression,
                        Name = InputList[[i]]$Name,
                        gName = "Cdk1", SpanVal = SpanVect[i], CatOrder = NULL)
  print(p)

}



for(i in 1:length(InputList)){

  p <- PlotOnPseudotime(WorkStruct = InputList[[i]]$OrderedData,
                        Expression = InputList[[i]]$Expression,
                        Name = InputList[[i]]$Name,
                        gName = "Cdk2", SpanVal = SpanVect[i], CatOrder = NULL)
  print(p)

}




for(i in 1:length(InputList)){

  p <- PlotOnPseudotime(WorkStruct = InputList[[i]]$OrderedData,
                        Expression = InputList[[i]]$Expression,
                        Name = InputList[[i]]$Name,
                        gName = "Cdk4", SpanVal = SpanVect[i], CatOrder = NULL)
  print(p)

}
















# Launch dashboard

# CompareAcrossData(InputList, NULL)






# AllFreman_Mouse
#
#
Genes.Buett <- Data.Buet$Analysis$Genes[[length(Data.Buet$Analysis$Genes)]]
Genes.Sasa <- Data.Sasa$Analysis$Genes[[length(Data.Sasa$Analysis$Genes)]]
Genes.Kowa <- Data.Kowa$Analysis$Genes[[length(Data.Kowa$Analysis$Genes)]]


Genes.Buett.Norm <- Data.Buet.Norm$Analysis$Genes[[length(Data.Buet.Norm$Analysis$Genes)]]
Genes.Sasa.Norm <- Data.Sasa.Norm$Analysis$Genes[[length(Data.Sasa.Norm$Analysis$Genes)]]
Genes.Kowa.Norm <- Data.Kowa.Norm$Analysis$Genes[[length(Data.Kowa.Norm$Analysis$Genes)]]

length(Genes.Buett)
length(Genes.Kowa)
length(Genes.Sasa)


length(Genes.Buett.Norm)
length(Genes.Kowa.Norm)
length(Genes.Sasa.Norm)



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

AllFreman_Mouse <- unique(AllFreman_Mouse)

AllFreman_Mouse_CCP <- c(Freeman_G1S_CC4_Mouse, Freeman_G1S_CC4A_Mouse, Freeman_G2_CC6A_Mouse,
                         Freeman_G2M_CC6_Mouse, Freeman_M_CC6B_Mouse, Freeman_S_CC4B_Mouse)


AllFreman_Mouse_CCP <- unique(AllFreman_Mouse_CCP)

AllWit_Mouse <- ConvertNames("human", "mouse", unlist(StageAssociation_Whit_Ext[-c(1,2)], use.names = FALSE))


WitAndFree <- intersect(AllFreman_Mouse, AllWit_Mouse)





library(VennDiagram)


#
# grid.newpage()
# draw.pairwise.venn(area1 = length(Genes.Buett), area2 = length(AllFreman_Mouse),
#                    cross.area = length(intersect(Genes.Buett, AllFreman_Mouse)),
#                    category = c("Buettner et al.", "Freeman"),
#                    fill = c("blue", "green"), cex = 3)
#
# grid.newpage()
# draw.pairwise.venn(area1 = length(Genes.Buett), area2 = length(AllWit_Mouse),
#                    cross.area = length(intersect(Genes.Buett, AllWit_Mouse)),
#                    category = c("Buettner et al.", "Withfield"),
#                    fill = c("blue", "green"), cex = 3)
#
#
# grid.newpage()
# draw.pairwise.venn(area1 = length(Genes.Buett), area2 = length(WitAndFree),
#                    cross.area = length(intersect(Genes.Buett, WitAndFree)),
#                    category = c("Buettner et al.", "Withfield"),
#                    fill = c("blue", "green"), cex = 3)
#
#
#
#
#
# grid.newpage()
# draw.quad.venn(area1 = length(Genes.Buett), area2 = length(Genes.Kowa),
#                area3 = length(Genes.Sasa), area4 = length(AllFreman_Mouse),
#                n12 = length(intersect(Genes.Buett, Genes.Kowa)),
#                n23 = length(intersect(Genes.Sasa, Genes.Kowa)),
#                n13 = length(intersect(Genes.Buett, Genes.Sasa)),
#                n14 = length(intersect(Genes.Buett, AllFreman_Mouse)),
#                n24 = length(intersect(AllFreman_Mouse, Genes.Kowa)),
#                n34 = length(intersect(Genes.Sasa, AllFreman_Mouse)),
#                n123 = length(intersect(Genes.Buett, intersect(Genes.Kowa, Genes.Sasa))),
#                n124 = length(intersect(Genes.Buett, intersect(Genes.Kowa, AllFreman_Mouse))),
#                n134 = length(intersect(Genes.Buett, intersect(Genes.Sasa, AllFreman_Mouse))),
#                n234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, AllFreman_Mouse))),
#                n1234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, intersect(Genes.Buett, AllFreman_Mouse)))),
#                category = c("Buettner et al.", "Kowalczyk et al.", "Sasagawa et al.", "Freeman"),
#                fill = c("red", "green", "blue", "white"), cex = 3)
# grid.newpage()
#
#
#
#
#
#
# grid.newpage()
# draw.quad.venn(area1 = length(Genes.Buett), area2 = length(Genes.Kowa),
#                area3 = length(Genes.Sasa), area4 = length(AllWit_Mouse),
#                n12 = length(intersect(Genes.Buett, Genes.Kowa)),
#                n23 = length(intersect(Genes.Sasa, Genes.Kowa)),
#                n13 = length(intersect(Genes.Buett, Genes.Sasa)),
#                n14 = length(intersect(Genes.Buett, AllWit_Mouse)),
#                n24 = length(intersect(AllWit_Mouse, Genes.Kowa)),
#                n34 = length(intersect(Genes.Sasa, AllWit_Mouse)),
#                n123 = length(intersect(Genes.Buett, intersect(Genes.Kowa, Genes.Sasa))),
#                n124 = length(intersect(Genes.Buett, intersect(Genes.Kowa, AllWit_Mouse))),
#                n134 = length(intersect(Genes.Buett, intersect(Genes.Sasa, AllWit_Mouse))),
#                n234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, AllWit_Mouse))),
#                n1234 = length(intersect(Genes.Sasa, intersect(Genes.Kowa, intersect(Genes.Buett, AllWit_Mouse)))),
#                category = c("Buettner et al.", "Kowalczyk et al.", "Sasagawa et al.", "Freeman"),
#                fill = c("red", "green", "blue", "white"), cex = 3)
# grid.newpage()
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





grid.newpage()
draw.triple.venn(area1 = length(Genes.Buett), area2 = length(Genes.Kowa),
               area3 = length(Genes.Sasa),
               n12 = length(intersect(Genes.Buett, Genes.Kowa)),
               n23 = length(intersect(Genes.Sasa, Genes.Kowa)),
               n13 = length(intersect(Genes.Buett, Genes.Sasa)),
               n123 = length(intersect(Genes.Buett, intersect(Genes.Kowa, Genes.Sasa))),
               category = c("Buettner et al.", "Kowalczyk et al.", "Sasagawa et al."),
               fill = c("red", "green", "blue"), cex = 3)


grid.newpage()
draw.triple.venn(area1 = length(Genes.Buett.Norm), area2 = length(Genes.Kowa.Norm),
                 area3 = length(Genes.Sasa.Norm),
                 n12 = length(intersect(Genes.Buett.Norm, Genes.Kowa.Norm)),
                 n23 = length(intersect(Genes.Sasa.Norm, Genes.Kowa.Norm)),
                 n13 = length(intersect(Genes.Buett.Norm, Genes.Sasa.Norm)),
                 n123 = length(intersect(Genes.Buett.Norm, intersect(Genes.Kowa.Norm, Genes.Sasa.Norm))),
                 category = c("Buettner et al. (Norm)", "Kowalczyk et al. (Norm)", "Sasagawa et al. (Norm)"),
                 fill = c("red", "green", "blue"), cex = 3)


grid.newpage()
draw.pairwise.venn(area1 = length(Genes.Buett),
                   area2 = length(Genes.Buett.Norm),
                   cross.area = length(intersect(Genes.Buett, Genes.Buett.Norm)),
                   category = c("Buettner et al.", "Buettner et al. (Norm)"),
                   fill = c("red", "blue"), cex = 2)


grid.newpage()
draw.pairwise.venn(area1 = length(Genes.Kowa),
                   area2 = length(Genes.Kowa.Norm),
                   cross.area = length(intersect(Genes.Kowa, Genes.Kowa.Norm)),
                   category = c("Kowalczyk et al.", "Kowalczyk et al. (Norm)"),
                   fill = c("red", "blue"), cex = 2)


grid.newpage()
draw.pairwise.venn(area1 = length(Genes.Sasa),
                   area2 = length(Genes.Sasa.Norm),
                   cross.area = length(intersect(Genes.Sasa, Genes.Sasa.Norm)),
                   category = c("Sasagawa et al.", "Sasagawa et al. (Norm)"),
                   fill = c("red", "blue"), cex = 2)





















AllSample <- intersect(Genes.Buett, intersect(Genes.Sasa, Genes.Kowa))
AllSample.Norm <- intersect(Genes.Buett.Norm, intersect(Genes.Sasa.Norm, Genes.Kowa.Norm))

grid.newpage()
draw.pairwise.venn(area1 = length(AllSample),
                   area2 = length(AllSample.Norm),
                   cross.area = length(intersect(AllSample, AllSample.Norm)),
                   category = c("Intersection", "Intersection (Norm)"),
                   fill = c("red", "blue"), cex = 2)



grid.newpage()
draw.triple.venn(area1 = length(AllSample), area2 = length(AllWit_Mouse),
                 area3 = length(AllFreman_Mouse),
                 n12 = length(intersect(AllSample, AllWit_Mouse)),
                 n23 = length(intersect(AllFreman_Mouse, AllWit_Mouse)),
                 n13 = length(intersect(AllSample, AllFreman_Mouse)),
                 n123 = length(intersect(AllSample, intersect(AllWit_Mouse, AllFreman_Mouse))),
                 category = c("Core module", "Withfield", "Giotti"),
                 fill = c("red", "green", "blue"), cex = 3)



grid.newpage()
draw.triple.venn(area1 = length(AllSample.Norm), area2 = length(AllWit_Mouse),
                 area3 = length(AllFreman_Mouse),
                 n12 = length(intersect(AllSample.Norm, AllWit_Mouse)),
                 n23 = length(intersect(AllFreman_Mouse, AllWit_Mouse)),
                 n13 = length(intersect(AllSample.Norm, AllFreman_Mouse)),
                 n123 = length(intersect(AllSample.Norm, intersect(AllWit_Mouse, AllFreman_Mouse))),
                 category = c("Core module", "Withfield", "Giotti"),
                 fill = c("red", "green", "blue"), cex = 3)


grid.newpage()

length(intersect(MouseGenes_GOCellCycle, AllFreman_Mouse)) / length(MouseGenes_GOCellCycle)
length(intersect(AllSample, AllFreman_Mouse)) / length(MouseGenes_GOCellCycle)



length(intersect(MouseGenes_GOCellCycle, AllWit_Mouse)) / length(MouseGenes_GOCellCycle)
length(intersect(AllSample, AllWit_Mouse)) / length(MouseGenes_GOCellCycle)







length(intersect(AllSample, MouseGenes_GOCellCycle))

New <- setdiff(AllSample, MouseGenes_GOCellCycle)


toupper(New) %in% Human_GOCellCycle









BaseBuett <- Genes.Buett
BaseKowa <- Genes.Kowa
BaseSasa <- Genes.Sasa


Buet.Sasa <- intersect(Genes.Buett, Genes.Sasa)
Buet.Kowa <- intersect(Genes.Buett, Genes.Kowa)
Sasa.Kowa <- intersect(Genes.Sasa, Genes.Kowa)

View(matrix(setdiff(unique(c(Buet.Sasa, Buet.Kowa, Sasa.Kowa)), AllSample), ncol = 1))



All.Buett <- length(AllSample)/length(BaseBuett)
Two.Buett <- length(setdiff(c(Buet.Sasa, Buet.Kowa), AllSample))/length(BaseBuett)
One.Buett <- length(setdiff(BaseBuett, c(BaseKowa, BaseSasa)))/length(BaseBuett)


All.Sasa <- length(AllSample)/length(BaseSasa)
Two.Sasa <- length(setdiff(c(Buet.Sasa, Sasa.Kowa), AllSample))/length(BaseSasa)
One.Sasa <- length(setdiff(BaseSasa, c(BaseBuett, BaseKowa)))/length(BaseSasa)


All.Kowa <- length(AllSample)/length(BaseKowa)
Two.Kowa <- length(setdiff(c(Sasa.Kowa, Buet.Kowa), AllSample))/length(BaseKowa)
One.Kowa <- length(setdiff(BaseKowa, c(BaseBuett, BaseSasa)))/length(BaseKowa)



barplot(c(One.Buett, Two.Buett, All.Buett), main="Buettner")
barplot(c(One.Kowa, Two.Kowa, All.Kowa), main="Kowa")
barplot(c(One.Sasa, Two.Sasa, All.Sasa), main="Sasa")


length(intersect(BaseSasa, MouseGenes_GOCellCycle))/length(BaseSasa)
length(intersect(BaseKowa, MouseGenes_GOCellCycle))/length(BaseKowa)
length(intersect(BaseBuett, MouseGenes_GOCellCycle))/length(BaseBuett)






# CatOrder = c("G0", "G0+G1", "G1", "G1+S", "S", "S+G2M", "G2M", "G2M+G0", "G1M+G1")



# DataStruct = Data.Buet
# ProcStruct = Proc.Exp.Buet
# MinNodes = 2
# FiltMax = .2
# Thr = .9
# QuantSel = .75
# SinglePeack = FALSE





View(matrix(sort(unique(AllSample)), ncol=8))









PeackGenes.Sasa <-  GetGenesWithPeaks(DataStruct = Data.Sasa,
                                      ProcStruct = Proc.Exp.Sasa,
                                      FiltMax = .1, Thr = .8,
                                      QuantSel = .5, SinglePeack = FALSE,
                                      AllGenes = TRUE)
PeackGenes.Kowa <-  GetGenesWithPeaks(DataStruct = Data.Kowa,
                                      ProcStruct = Proc.Exp.Kowa,
                                      FiltMax = .1, Thr = .8,
                                      QuantSel = .5, SinglePeack = FALSE,
                                      AllGenes = TRUE)
PeackGenes.Buet <-  GetGenesWithPeaks(DataStruct = Data.Buet,
                                      ProcStruct = Proc.Exp.Buet,
                                      FiltMax = .1, Thr = .8,
                                      QuantSel = .5, SinglePeack = FALSE,
                                      AllGenes = TRUE)










PeackGenes.Sasa.Norm <-  GetGenesWithPeaks(DataStruct = Data.Sasa.Norm,
                                      ProcStruct = Proc.Exp.Sasa.Norm,
                                      FiltMax = .1, Thr = .8,
                                      QuantSel = .5, SinglePeack = FALSE,
                                      AllGenes = TRUE)
PeackGenes.Kowa.Norm <-  GetGenesWithPeaks(DataStruct = Data.Kowa.Norm,
                                      ProcStruct = Proc.Exp.Kowa.Norm,
                                      FiltMax = .1, Thr = .8,
                                      QuantSel = .5, SinglePeack = FALSE,
                                      AllGenes = TRUE)
PeackGenes.Buet.Norm <-  GetGenesWithPeaks(DataStruct = Data.Buet.Norm,
                                      ProcStruct = Proc.Exp.Buet.Norm,
                                      FiltMax = .1, Thr = .8,
                                      QuantSel = .5, SinglePeack = FALSE,
                                      AllGenes = TRUE)



# DataStruct = Data.Buet
# ProcStruct = Proc.Exp.Buet
# FiltMax = .1
# Thr = .9
# QuantSel = .5
# SinglePeack = FALSE
# MinNodes = 2










names(PeackGenes.Kowa)
names(PeackGenes.Buet)
names(PeackGenes.Sasa)


G1.kowa <- unlist(PeackGenes.Kowa[grep("G1", names(PeackGenes.Kowa))])
G1.buet <- unlist(PeackGenes.Buet[grep("G1", names(PeackGenes.Buet))])
G1.sasa <- unlist(PeackGenes.Sasa[grep("G1", names(PeackGenes.Sasa))])

S.kowa <- unlist(PeackGenes.Kowa[grep("S", names(PeackGenes.Kowa))])
S.buet <- unlist(PeackGenes.Buet[grep("S", names(PeackGenes.Buet))])
S.sasa <- unlist(PeackGenes.Sasa[grep("S", names(PeackGenes.Sasa))])

G2M.kowa <- unlist(PeackGenes.Kowa[grep("G2/M", names(PeackGenes.Kowa))])
G2M.buet <- unlist(PeackGenes.Buet[grep("G2M", names(PeackGenes.Buet))])
G2M.sasa <- unlist(PeackGenes.Sasa[grep("G2/M", names(PeackGenes.Sasa))])


G1.kowa.Norm <- unlist(PeackGenes.Kowa.Norm[grep("G1", names(PeackGenes.Kowa.Norm))])
G1.buet.Norm <- unlist(PeackGenes.Buet.Norm[grep("G1", names(PeackGenes.Buet.Norm))])
G1.sasa.Norm <- unlist(PeackGenes.Sasa.Norm[grep("G1", names(PeackGenes.Sasa.Norm))])

S.kowa.Norm <- unlist(PeackGenes.Kowa.Norm[grep("S", names(PeackGenes.Kowa.Norm))])
S.buet.Norm <- unlist(PeackGenes.Buet.Norm[grep("S", names(PeackGenes.Buet.Norm))])
S.sasa.Norm <- unlist(PeackGenes.Sasa.Norm[grep("S", names(PeackGenes.Sasa.Norm))])

G2M.kowa.Norm <- unlist(PeackGenes.Kowa.Norm[grep("G2/M", names(PeackGenes.Kowa.Norm))])
G2M.buet.Norm <- unlist(PeackGenes.Buet.Norm[grep("G2M", names(PeackGenes.Buet.Norm))])
G2M.sasa.Norm <- unlist(PeackGenes.Sasa.Norm[grep("G2/M", names(PeackGenes.Sasa.Norm))])








# G0.kowa <- unlist(PeackGenes.Kowa[c(1,10)])

# G1.kowa <- setdiff(unlist(PeackGenes.Kowa[2:6]),
#                    unlist(PeackGenes.Kowa[-c(2:6)]))
# G1.buet <- setdiff(unlist(PeackGenes.Buet[1:2]),
#                    unlist(PeackGenes.Buet[-c(1:2)]))
# G1.sasa <- setdiff(unlist(PeackGenes.Sasa[1:2]),
#                    unlist(PeackGenes.Sasa[-c(1:2)]))
#
# S.kowa <- setdiff(unlist(PeackGenes.Kowa[6:8]),
#                   unlist(PeackGenes.Kowa[-c(6:8)]))
# S.buet <- setdiff(unlist(PeackGenes.Buet[2:4]),
#                   unlist(PeackGenes.Buet[-c(2:4)]))
# S.sasa <- setdiff(unlist(PeackGenes.Sasa[2:4]),
#                   unlist(PeackGenes.Sasa[-c(2:4)]))
#
# G2M.kowa <- setdiff(unlist(PeackGenes.Kowa[8:10]),
#                     unlist(PeackGenes.Kowa[-c(8:10)]))
# G2M.buet <- setdiff(unlist(PeackGenes.Buet[4:6]),
#                     unlist(PeackGenes.Buet[-c(4:6)]))
# G2M.sasa <- setdiff(unlist(PeackGenes.Sasa[4:6]),
#                     unlist(PeackGenes.Sasa[-c(4:6)]))

# setdiff(G1.kowa, S.kowa)

# setdiff(PeackGenes.Buet[[2]], PeackGenes.Buet[[2]])



# G0.All <- G0.kowa

G1.All <- unique(c(G1.kowa, G1.buet, G1.sasa))

S.All <- c(S.kowa, S.buet, S.sasa)

G2M.All <- c(G2M.kowa, G2M.buet, G2M.sasa)

length(G1.All)
length(S.All)
length(G2M.All)






G1.Combined <- rbind(G1.All %in% G1.buet,
                     G1.All %in% G1.kowa,
                     G1.All %in% G1.sasa)

colnames(G1.Combined) <- G1.All
rownames(G1.Combined) <- c("Buet", "Kowa", "Sasa")

table(colSums(G1.Combined))
table(colSums(G1.Combined[-2,]))

table(colSums(G1.Combined[-1,]))
table(colSums(G1.Combined[-3,]))







grid.newpage()
draw.triple.venn(area1 = length(G1.buet), area2 = length(G1.sasa),
                 area3 = length(G1.kowa), scaled = TRUE,
                 n12 = length(intersect(G1.buet, G1.sasa)),
                 n23 = length(intersect(G1.kowa, G1.sasa)),
                 n13 = length(intersect(G1.buet, G1.kowa)),
                 n123 = length(intersect(G1.buet, intersect(G1.sasa, G1.kowa))),
                 category = c("Buettner et al. (G1)", "Sasagawa et al. (G1)", "Kowalczyk et al. (G1)"),
                 fill = c("red", "green", "blue"), cex = 3, main = "pippo")




grid.newpage()
draw.triple.venn(area1 = length(S.buet), area2 = length(S.sasa),
                 area3 = length(S.kowa), scaled = TRUE,
                 n12 = length(intersect(S.buet, S.sasa)),
                 n23 = length(intersect(S.kowa, S.sasa)),
                 n13 = length(intersect(S.buet, S.kowa)),
                 n123 = length(intersect(S.buet, intersect(S.sasa, S.kowa))),
                 category = c("Buettner et al. (S)", "Sasagawa et al. (S)", "Kowalczyk et al. (S)"),
                 fill = c("red", "green", "blue"), cex = 3, main = "pippo")



grid.newpage()
draw.triple.venn(area1 = length(G2M.buet), area2 = length(G2M.sasa),
                 area3 = length(G2M.kowa), scaled = TRUE,
                 n12 = length(intersect(G2M.buet, G2M.sasa)),
                 n23 = length(intersect(G2M.kowa, G2M.sasa)),
                 n13 = length(intersect(G2M.buet, G2M.kowa)),
                 n123 = length(intersect(G2M.buet, intersect(G2M.sasa, G2M.kowa))),
                 category = c("Buettner et al. (G2M)", "Sasagawa et al. (G2M)", "Kowalczyk et al. (G2M)"),
                 fill = c("red", "green", "blue"), cex = 3, main = "pippo")








grid.newpage()
draw.triple.venn(area1 = length(G1.buet.Norm), area2 = length(G1.sasa.Norm),
                 area3 = length(G1.kowa.Norm), scaled = TRUE,
                 n12 = length(intersect(G1.buet.Norm, G1.sasa.Norm)),
                 n23 = length(intersect(G1.kowa.Norm, G1.sasa.Norm)),
                 n13 = length(intersect(G1.buet.Norm, G1.kowa.Norm)),
                 n123 = length(intersect(G1.buet.Norm, intersect(G1.sasa.Norm, G1.kowa.Norm))),
                 category = c("Buettner et al. (Norm - G1)", "Sasagawa et al. (Norm - G1)", "Kowalczyk et al. (Norm - G1)"),
                 fill = c("red", "green", "blue"), cex = 3, main = "pippo")




grid.newpage()
draw.triple.venn(area1 = length(S.buet.Norm), area2 = length(S.sasa.Norm),
                 area3 = length(S.kowa.Norm), scaled = TRUE,
                 n12 = length(intersect(S.buet.Norm, S.sasa.Norm)),
                 n23 = length(intersect(S.kowa.Norm, S.sasa.Norm)),
                 n13 = length(intersect(S.buet.Norm, S.kowa.Norm)),
                 n123 = length(intersect(S.buet.Norm, intersect(S.sasa.Norm, S.kowa.Norm))),
                 category = c("Buettner et al. (S)", "Sasagawa et al. (S)", "Kowalczyk et al. (S)"),
                 fill = c("red", "green", "blue"), cex = 3, main = "pippo")



grid.newpage()
draw.triple.venn(area1 = length(G2M.buet.Norm), area2 = length(G2M.sasa.Norm),
                 area3 = length(G2M.kowa.Norm), scaled = TRUE,
                 n12 = length(intersect(G2M.buet.Norm, G2M.sasa.Norm)),
                 n23 = length(intersect(G2M.kowa.Norm, G2M.sasa.Norm)),
                 n13 = length(intersect(G2M.buet.Norm, G2M.kowa.Norm)),
                 n123 = length(intersect(G2M.buet.Norm, intersect(G2M.sasa.Norm, G2M.kowa.Norm))),
                 category = c("Buettner et al. (G2M)", "Sasagawa et al. (G2M)", "Kowalczyk et al. (G2M)"),
                 fill = c("red", "green", "blue"), cex = 3, main = "pippo")






sum(G2M.sasa %in% G2M.kowa)

G2M_Strong <- intersect(intersect(G2M.sasa, G2M.buet), G2M.kowa)
G1_Strong <- intersect(intersect(G1.sasa, G1.buet), G1.kowa)
S_Strong <- intersect(intersect(S.sasa, S.buet), S.kowa)

G2M_Strong.Norm <- intersect(intersect(G2M.sasa.Norm, G2M.buet.Norm), G2M.kowa.Norm)
G1_Strong.Norm <- intersect(intersect(G1.sasa.Norm, G1.buet.Norm), G1.kowa.Norm)
S_Strong.Norm <- intersect(intersect(S.sasa.Norm, S.buet.Norm), S.kowa.Norm)




fileConn <- file("~/Desktop/PriceCyc_Mouse.gmt")
writeLines(c(paste(c("G1", "G1 Strong", G1_Strong), collapse ='\t'),
             paste(c("S", "S Strong", S_Strong), collapse ='\t'),
             paste(c("G2M", "G2M Strong", G2M_Strong), collapse ='\t'),
             paste(c("Core", "Core module", AllSample), collapse ='\t')),
           fileConn)
close(fileConn)




fileConn <- file("~/Desktop/PriceCyc_Hum.gmt")
writeLines(c(paste(c("G1", "G1 Strong", toupper(G1_Strong)), collapse ='\t'),
             paste(c("S", "S Strong", toupper(S_Strong)), collapse ='\t'),
             paste(c("G2M", "G2M Strong", toupper(G2M_Strong)), collapse ='\t'),
             paste(c("Core", "Core module", toupper(AllSample)), collapse ='\t')),
           fileConn)
close(fileConn)








# G1.Sel <- rbind(G1.All %in% G1.kowa,
#                 G1.All %in% G1.buet,
#                 G1.All %in% G1.sasa)
# G1.All <- G1.All[which(colSums(G1.Sel) > 1)]
#
# S.Sel <- rbind(S.All %in% S.kowa,
#                 S.All %in% S.buet,
#                 S.All %in% S.sasa)
# S.All <- S.All[which(colSums(S.Sel) > 1)]
#
# G2M.Sel <- rbind(G2M.All %in% G2M.kowa,
#                 G2M.All %in% G2M.buet,
#                 G2M.All %in% G2M.sasa)
# G2M.All <- G2M.All[which(colSums(G2M.Sel) > 1)]

length(G1.All)
length(S.All)
length(G2M.All)


ToRemove <- unique(c(intersect(G1.All, S.All),
                     intersect(S.All, G2M.All),
                     intersect(G1.All, G2M.All))
)


# G1.Sel <- G1.All

G1.Sel <- G1.All
G1.Sel <- G1.Sel[!(G1.Sel %in% ToRemove)]

S.Sel <- S.All
S.Sel <- S.Sel[!(S.Sel %in% ToRemove)]

G2M.Sel <- G2M.All
G2M.Sel <- G2M.Sel[!(G2M.Sel %in% ToRemove)]

# G0.Sel <- G0.All
# G0.Sel <- G0.Sel[!(G0.Sel %in% ToRemove)]



sum(G1.buet %in% G1.Sel)
sum(S.buet %in% S.Sel)
sum(G2M.buet %in% G2M.Sel)

sum(G1.sasa %in% G1.Sel)
sum(S.sasa %in% S.Sel)
sum(G2M.sasa %in% G2M.Sel)

sum(G1.kowa %in% G1.Sel)
sum(S.kowa %in% S.Sel)
sum(G2M.kowa %in% G2M.Sel)





sum(G1.buet %in% G1.All)
sum(S.buet %in% S.All)
sum(G2M.buet %in% G2M.All)

sum(G1.sasa %in% G1.All)
sum(S.sasa %in% S.All)
sum(G2M.sasa %in% G2M.All)

sum(G1.kowa %in% G1.All)
sum(S.kowa %in% S.All)
sum(G2M.kowa %in% G2M.All)








DF <- data.frame(lapply(list(G1 = G1_Strong, S = S_Strong, G2M = G2M_Strong), length))

library(ggplot2)
library(reshape)

melt(DF)

ggplot(data = melt(DF), mapping = aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = 'identity') + guides(fill='none') +
  labs(title="Biomarker genes", y="Number of genes", x="Cell cycle phase")


# SelStageInfo <- list(G1 = G1.All, S = S.All, G2M = G2M.All)
# SelStageInfo <- list(G1 = G1.Sel, S = S.Sel, G2M = G2M.Sel)
# SelStageInfo <- list(G1 = G1_Strong, S = S_Strong, G2M = G2M_Strong)
SelStageInfo <- list(G1 = G1_Strong.Norm, S = S_Strong.Norm, G2M = G2M_Strong.Norm)

# SelStageInfo <- list(G1 = G1_Strong.Norm[!(G1_Strong.Norm %in% c(S_Strong.Norm, G2M_Strong.Norm))],
#                      S = S_Strong.Norm[!(S_Strong.Norm %in% c(G1_Strong.Norm, G2M_Strong.Norm))],
#                      G2M = G2M_Strong.Norm[!(G2M_Strong.Norm %in% c(G1_Strong.Norm, S_Strong.Norm))])

Staged.Sasa <- StageWithPeaks(DataStruct = Data.Sasa,
               ProcStruct = Proc.Exp.Sasa,
               ComputeG0 = FALSE,
               FiltMax = 0,
               Thr = .8,
               QuantSel = .8,
               StageInfo = SelStageInfo,
               MinNodes = 2,
               Mode = 1,
               G0Level = .9,
               Title = 'Sasagawa et al')

#
# TB <- table(Staged.Sasa$CellStages, Data.Sasa$Cats)
# t(TB)/colSums(TB)
#
# Data.Sasa$Analysis
#





Staged.Buet <- StageWithPeaks(DataStruct = Data.Buet,
                              ProcStruct = Proc.Exp.Buet,
                              FiltMax = .1,
                              Thr = .8,
                              QuantSel = .8,
                              StageInfo = SelStageInfo,
                              ComputeG0 = FALSE,
                              MinNodes = 2,
                              Mode = 1,
                              G0Level = .9,
                              Title = 'Buettner et al')









Staged.Kowa <- StageWithPeaks(DataStruct = Data.Kowa,
                              ProcStruct = Proc.Exp.Kowa,
                              FiltMax = .1,
                              Thr = .8,
                              QuantSel = .8,
                              StageInfo = SelStageInfo,
                              ComputeG0 = TRUE,
                              MinNodes = 2,
                              Mode = 1,
                              G0Level = 1.1,
                              Title = "Kowalczyk et al")


#
#
# library(GSEABase)
#
#
#
#
# devtools::install_url("https://github.com/satijalab/seurat/releases/download/v1.4.0/Seurat_1.4.0.16.tgz", binary = TRUE)
# library(Seurat)
#
#


#
#
# DataStruct = Data.Sasa
# ProcStruct = Proc.Exp.Sasa
# FiltMax = .1
# Thr = .75
# QuantSel = .75
# StageInfo = SelStageInfo
# ComputeG0 = TRUE
# MinNodes = 2
# Mode = 3
# G0Level = .7







#
# Staged.Kowa <- StageWithPeaks(DataStruct = Data.Kowa,
#                               ProcStruct = Proc.Exp.Kowa,
#                               FiltMax = .1,
#                               Thr = .75,
#                               QuantSel = .75,
#                               SinglePeack = TRUE,
#                               StageInfo = SelStageInfo,
#                               ComputeG0 = TRUE,
#                               MinNodes = 2,
#                               Mode = 3,
#                               G0Level = 1)
#







length(G0.Sel)
length(G1.Sel)
length(S.Sel)
length(G2M.Sel)


table(colSums(rbind(G1.All %in% G1.buet,
              G1.All %in% G1.sasa,
              G1.All %in% G1.kowa)))




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
