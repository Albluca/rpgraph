

# Check directionality --------------------------------------------


CheckDirectionality <- function(Data, Results, OrderClass, Depth = NULL){

  TaxonList <- getTaxonMap(Graph = makeGraph(Results), Data = Data)

  if(min(Results$Edges) == 0){
    Results$Edges <- Results$Edges + 1
  }

  Net <- graph.empty(n = length(unique(as.vector(Results$Edges))), directed = FALSE)
  V(Net)$name <- paste("V_", unique(as.vector(Results$Edges)), sep = '')

  for (i in 1:nrow(Results$Edges)) {
    Net <- add.edges(graph = Net, paste("V_", Results$Edges[i,], sep = ''))
  }

  EdgDir <- NULL

  if(is.null(Depth)){
    Depth <- vcount(Net)
  }

  for (Edg in E(Net)) {

    tNet <- delete_edges(Net, Edg)

    Verts <- V(Net)$name[get.edges(Net, Edg)]

    Nei <- neighborhood(graph = tNet, order = Depth, nodes = Verts)
    VertSeq1 <- unlist(strsplit(V(tNet)$name[Nei[[1]]], split = "V_"))
    VertSeq1 <- as.numeric(VertSeq1[-seq(from = 1, by = 2, along.with = VertSeq1)])

    VertSeq2 <- unlist(strsplit(V(tNet)$name[Nei[[2]]], split = "V_"))
    VertSeq2 <- as.numeric(VertSeq2[-seq(from = 1, by = 2, along.with = VertSeq2)])

    Clean_Levels <- unlist(TaxonList[VertSeq1])
    Order1 <- OrderClass[Clean_Levels[!is.na(Clean_Levels)]]

    Clean_Levels <- unlist(TaxonList[VertSeq2])
    Order2 <- OrderClass[Clean_Levels[!is.na(Clean_Levels)]]

    if(sum(is.finite(Order1)) < 3 | sum(is.finite(Order2)) < 3 ){
      EdgDir <- rbind(EdgDir, c(Verts, NA, NA))
      next()
    }

    WT <- wilcox.test(x = Order1, y = Order2, na.rm=TRUE)

    if(median(Order1, na.rm = TRUE) < median(Order2, na.rm = TRUE)){
      EdgDir <- rbind(EdgDir, c(Verts, 1, WT$p.value))
      next()
    }

    if(median(Order1, na.rm = TRUE) > median(Order2, na.rm = TRUE)){
      EdgDir <- rbind(EdgDir, c(Verts, 2, WT$p.value))
      next()
    }

    if(median(Order1, na.rm = TRUE) == median(Order2, na.rm = TRUE)){
      EdgDir <- rbind(EdgDir, c(Verts, 0, WT$p.value))
      next()
    }

  }

  colnames(EdgDir) <- c("Source", "Target", "Direction", "P.val")

  return(data.frame(EdgDir, stringsAsFactors = FALSE))

}


GetDirectionalityIndex <- function(InfoData) {

  VertDist <- distances(InfoData$Net, mode = 'out')
  Dists <- cbind(which(lower.tri(VertDist), arr.ind = TRUE), VertDist[lower.tri(VertDist)], t(VertDist)[lower.tri(t(VertDist))])
  Dists <- cbind(V(InfoData$Net)$name[Dists[,1]], V(InfoData$Net)$name[Dists[,2]], Dists[,3:4])
  colnames(Dists) <- c("V1", "V2", "V1->V2", "V2->V1")

  DirectionalityIndex <- (sum(is.infinite(as.numeric(Dists[,3])) &
                                !is.infinite(as.numeric(Dists[,4]))) +
                            sum(!is.infinite(as.numeric(Dists[,3])) &
                                  is.infinite(as.numeric(Dists[,4]))))/sum(!is.infinite(as.numeric(Dists[,3])) |
                                                                             !is.infinite(as.numeric(Dists[,4])))

  return(DirectionalityIndex)
}






# CompareNets -------------------------------------------------------------


CompareNet <- function(G1, G2, RemNodes = 2, Tries = 10000, DoIso = FALSE) {

  if(DoIso){
    Full_Iso <- graph.get.isomorphisms.vf2(graph1 = G1, graph2 = G2)

    if(length(Full_Iso)>0){
      return(0)
    }
  }


  pb <- txtProgressBar(min = 1, max = Tries, initial = 1, style = 3)

  for (Retries in 1:Tries) {
    setTxtProgressBar(pb, Retries)
    RemVert <- sample(x = V(G1), size = RemNodes, replace = FALSE)
    tNet <- delete_vertices(G1, RemVert)
    Part_Iso <- graph.get.subisomorphisms.vf2(graph1 = G2, graph2 = tNet)
    if(length(Part_Iso)>0){
      RetVal <- list(rem)
      close(pb)
      return(list(SubIso = Part_Iso, RemVert = RemVert))
    }
  }

  close(pb)
  return(NULL)

}





