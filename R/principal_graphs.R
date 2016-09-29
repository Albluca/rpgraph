# File description --------------------------------------------
#
# This fine contains functions designed to interaqct with VDAOEngine.jar
#
# Only the main function will be exported and made available to the user
#
#



# Configuration Functions --------------------------------------------


#' Define the grammar to produde a circular principal graph
#'
#' This is an internal function and should not be accessed directly by the user.
#'
#' @return A parameter set that can be used to construct a circular principal graph
CircleConfiguration <- function(){

  params <- NULL

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a closed reactangle spanned by first two PCs
  params$initalgorithm = 2

  # universal parameters (should be rarely modified)
  params$epochs[[1]]$id = 1
  params$epochs[[1]]$minimize = 1
  params$epochs[[1]]$eps = 0.01
  params$epochs[[1]]$epsSLAU = 0.01
  params$epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  params$epochs[[1]]$grammartype = 'circle';

  # elasticity parameters
  params$epochs[[1]]$ep = 0.01
  params$epochs[[1]]$rp = 0.1

  # default number of nodes
  params$epochs[[1]]$numiter = 20

  return(params)
}

#' Define the grammar to produde a tree principal graph
#'
#' This is an internal function and should not be accessed directly by the user.
#'
#' @return A parameter set that can be used to construct a tree principal graph
DefaultPrincipalTreeConfiguration <- function(){

  params <- NULL

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a piece of first principal component
  params$initalgorithm = 0

  # universal parameters (should be rarely modified)
  params$epochs[[1]]$id = 1
  params$epochs[[1]]$minimize = 1
  params$epochs[[1]]$eps = 0.01
  params$epochs[[1]]$epsSLAU = 0.01
  params$epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  params$epochs[[1]]$grammartype = 'treeWithPruning';

  # elasticity parameters
  params$epochs[[1]]$ep = 0.01
  params$epochs[[1]]$rp = 0.1

  # default number of nodes
  params$epochs[[1]]$numiter = 20

  return(params)


}

#' Define the grammar to produde a curve (a path) principal graph
#'
#' This is an internal function and should not be accessed directly by the user.
#'
#' @return A parameter set that can be used to construct a curve (a path) principal graph
CurveConfiguration <- function(){

  params <- NULL

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a piece of first principal component
  params$initalgorithm = 0

  # universal parameters (should be rarely modified)
  params$epochs[[1]]$id = 1
  params$epochs[[1]]$minimize = 1
  params$epochs[[1]]$eps = 0.01
  params$epochs[[1]]$epsSLAU = 0.01
  params$epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  params$epochs[[1]]$grammartype = 'curve'

  # elasticity parameters
  params$epochs[[1]]$ep = 0.01
  params$epochs[[1]]$rp = 0.1

  # default number of nodes
  params$epochs[[1]]$numiter = 20

  return(params)
}

#' Define the grammar to produde a robust tree principal graph
#'
#' This is an internal function and should not be accessed directly by the user
#'
#' @return A parameter set that can be used to construct a robust tree principal graph
RobustPrincipalTreeConfiguration <- function(){

  params <- NULL

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a piece of first principal component
  params$initalgorithm = 1

  # universal parameters (should be rarely modified)
  params$epochs[[1]]$id = 1
  params$epochs[[1]]$minimize = 1
  params$epochs[[1]]$eps = 0.01
  params$epochs[[1]]$epsSLAU = 0.01
  params$epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  params$epochs[[1]]$grammartype = 'treeWithPruning';

  # elasticity parameters
  params$epochs[[1]]$ep = 0.01
  params$epochs[[1]]$rp = 0.1

  # default number of nodes
  params$epochs[[1]]$numiter = 30

  # robustness parameters
  params$epochs[[1]]$robust = TRUE
  params$epochs[[1]]$trimradius = 0.1

  return(params)
}

#' Define the grammar to produde a tree principal graph using a two stage approach.
#'
#' This is an internal function and should not be accessed directly by the user.
#' This function is still experimental and not is not working yet.
#'
#' @return A parameter set that can be used to construct a robust tree principal graph
TwoStagesRobustPrincipalTreeConfiguration <- function(){

  params <- list()

  params$algtype = 'grammar'
  params$initstretch1 = 1

  # initialization by a piece of first principal component
  params$initalgorithm = 0

  ###########
  # Epoch 1 #
  ###########

  epochs <- list()
  epochs[[1]] <- NA

  # universal parameters (should be rarely modified)
  epochs[[1]]$id = 1
  epochs[[1]]$minimize = 1
  epochs[[1]]$eps = 0.01
  epochs[[1]]$epsSLAU = 0.01
  epochs[[1]]$numiterSLAU = 1000

  # what grammar to use
  epochs[[1]]$grammartype = 'treeWithPruning';

  # elasticity parameters
  epochs[[1]]$ep = 0.01
  epochs[[1]]$rp = 0.1

  # default number of nodes
  epochs[[1]]$numiter = 30

  # robustness parameters
  epochs[[1]]$robust = FALSE
  epochs[[1]]$trimradius = 0.1

  ###########
  # Epoch 2 #
  ###########

  epochs[[2]] <- NA

  # universal parameters (should be rarely modified)
  epochs[[2]]$id = 1
  epochs[[2]]$minimize = 1
  epochs[[2]]$eps = 0.01
  epochs[[2]]$epsSLAU = 0.01
  epochs[[2]]$numiterSLAU = 1000

  # what grammar to use
  epochs[[2]]$grammartype = 'treeWithPruning';

  # elasticity parameters
  epochs[[2]]$ep = 0.0001
  epochs[[2]]$rp = 0.0001

  # default number of nodes
  epochs[[2]]$numiter = 30

  # robustness parameters
  epochs[[2]]$robust = TRUE
  epochs[[2]]$trimradius = 0.1

  params$epochs <- epochs

  return(params)
}










# Data extraction Functions --------------------------------------------


#' Produce a graph structure from the coputed elastic principal graph details
#'
#' The function depends on a working JVM with the VDAOEngine.jar library in the libpath.
#' This is an internal function and should not be accessed directly by the user
#'
#' @param PrintGraph An elastic principal graph structure as returned from \code{link{computeElasticPrincipalGraph}}
#'
#' @return A graph structure
makeGraph <- function(PrintGraph){

  Graph = .jnew('vdaoengine/analysis/grammars/Graph')

  Utils <- .jnew(class = 'vdaoengine/analysis/grammars/Utils')

  Utils$setNodesForGraph(Graph, .jarray(.jfloat(data.matrix(PrintGraph$Nodes)),dispatch=T))
  if(min(PrintGraph$Edges) == 0){
    Utils$setEdgesForGraph(Graph, .jarray(matrix(as.integer(PrintGraph$Edges), ncol = 2),dispatch=T))
  } else {
    Utils$setEdgesForGraph(Graph, .jarray(matrix(as.integer(PrintGraph$Edges-1), ncol = 2),dispatch=T))
  }

  Graph$defineStarsFromPrimitiveGraphStructure()
  Graph$compileNodesInArrays()

  return(Graph)

}

#' Project a set of datapoints onto the nodes of principal graph
#'
#' The function depends on a working JVM with the VDAOEngine.jar library in the libpath.
#'
#' @param PrintGraph An elastic principal graph structure as returned from \code{link{computeElasticPrincipalGraph}}
#' @param Data A matrix containing a set of points with the correct number of dimensions. Each row represent a point
#'
#' @return A list. Each elements of the list represent the vertex of the principal graph and contains a vector reporting
#' the cells (as row number) that are associted with that particular vertex. The vector equal to NA indicates that no
#' cells are associated with that vertex;
#' @export
getTaxonMap <- function(Results, Data){

  Graph <- makeGraph(Results)
  
  NumberOfNodes <- Graph$Nodes$size()

  numpoints <- as.integer(nrow(Data))
  coordnum <- as.integer(ncol(Data))

  Dataset <- .jnew('vdaoengine/data/VDataSet')
  Dataset$massif <- .jarray(.jfloat(data.matrix(Data)),dispatch=T)
  Dataset$pointCount = numpoints
  Dataset$coordCount = coordnum

  Elo = .jnew(class = 'vdaoengine/analysis/grammars/ElasticEnergyOptimization', Dataset, Graph)

  Elo$calcTaxons()

  TaxonMap <- list()

  for (i in 1:NumberOfNodes) {
    Tx <- Elo$taxons$get(as.integer(i-1))
    TaxonSize <- Tx$size()

    if(TaxonSize == 0){
      TaxonVect <- NA
    } else {
      TaxonVect <-NULL
      for (j in 1:TaxonSize) {
        TaxonVect[j] <- Tx$get(as.integer(j-1))+1
      }
    }

    TaxonMap[[i]] <- TaxonVect

  }

  return(TaxonMap)

}






# computeElPT --------------------------------------------





#' Compute the the elastic principal graph
#'
#' The function depends on a working JVM with the VDAOEngine.jar library in the libpath.
#' This is an internal function and should not be accessed directly by the user.
#'
#' @param Data The data points to be used to compute the graph
#' @param NumNodes The number of vertices to be used to construct the graph
#' @param Parameters The name of a parameter function
#' @param ... Additional parameter to be procesed by the Java procedure. Currently only EP, RP and TrimRadius are supported
#'
#' @return A list describing the principal elastic graph
computeElPT <- function(Data, NumNodes, Parameters, ...){

  StartTime <- Sys.time()

  EP <- NULL
  RP <- NULL
  TrimRadius <- -1

  AddPar <- list(...)

  if(length(AddPar)>0){

    print("Additional parameters found")
    print(AddPar)

  }

  if(length(AddPar)>0){
    if(any(names(AddPar)=="EP")){
      EP <- AddPar[["EP"]]
      print("EP set from user")
    }

    if(any(names(AddPar)=="RP")){
      RP <- AddPar[["RP"]]
      print("RP set from user")
    }

    if(any(names(AddPar)=="TrimRadius")){
      TrimRadius <- AddPar[["TrimRadius"]]
      print("TrimRadius set from user")
    }

  }


  cat("Configuring engine .")

  Config <- .jnew(class = "vdaoengine/analysis/grammars/ConfigFile")

  Config$algtype <- Parameters$algtype

  tArray <- Config$stretchInitCoeffs
  tArray[1] <- Parameters$initstretch1

  Config$stretchInitCoeffs <- .jfloat(tArray)

  Config$initStrategy <- as.integer(Parameters$initalgorithm)

  cat(".")

  for (i in 1:length(Parameters$epochs)) {

    Epoch <- .jnew(class = "vdaoengine/analysis/elmap/ElmapAlgorithmEpoch")

    Epoch$grammarType <- Parameters$epochs[[i]]$grammartype

    Epoch$EP <- .jfloat(Parameters$epochs[[i]]$ep)

    if(!is.null(EP)){
      if(EP[i]>0){
        Epoch$EP <- .jfloat(EP[i])
      }
    }

    Epoch$RP = .jfloat(Parameters$epochs[[i]]$rp)

    if(!is.null(RP)){
      if(RP[i]>0){
        Epoch$RP <- .jfloat(RP[i])
      }
    }

    Epoch$numberOfIterations <- as.integer(NumNodes[i])
    Epoch$maxNumberOfIterationsForSLAU <- as.integer(Parameters$epochs[[i]]$numiterSLAU)
    Epoch$epsConvergence <- .jfloat(Parameters$epochs[[i]]$eps)
    Epoch$epsConvergenceSLAU <- .jfloat(Parameters$epochs[[i]]$epsSLAU)

    if(!is.null(Parameters$epochs[[1]]$robust)){
      Epoch$robust <- Parameters$epochs[[i]]$robust
    }

    if(!is.null(Parameters$epochs[[1]]$trimradius)){
      Epoch$trimradius <- .jfloat(Parameters$epochs[[i]]$trimradius)
    }

    if(TrimRadius>0){
      Epoch$robust = TRUE
      Epoch$trimradius = TrimRadius[i]
    }

    Config$epochs$add(Epoch)

  }

  cat(".")

  cpg <- .jnew(class = "vdaoengine/analysis/grammars/ComputePrincipalGraph")

  cpg$config <- Config

  cat(".")

  cpg$setDataSetAsMassif(.jarray(.jfloat(data.matrix(Data)),dispatch=T))
  cat(".")

  print("")
  print("Running engine")

  report <- read.delim(text = cpg$compute())

  NodeSize <- cpg$graph$countNumberOfPointsProjected(.jcast(cpg$dataset, new.class = "/java/lang/String" ))

  NodePositions <- cpg$graph$getNodePositions()
  NodePositions <- .jevalArray(NodePositions, simplify = TRUE)

  Edges <- cpg$graph$getEdgeTable()
  Edges <- .jevalArray(Edges, simplify = TRUE)

  EndTime <- Sys.time()


  StructuredReturn <- list()
  StructuredReturn$Nodes <- NodePositions
  StructuredReturn$Edges <- Edges
  StructuredReturn$Report <- report
  StructuredReturn$EndTime <- EndTime
  StructuredReturn$StartTime <- StartTime
  StructuredReturn$NodeSize <- NodeSize


  return(StructuredReturn)
}




# computeElasticPrincipalGraph --------------------------------------------


#' Compute the elastic principal graph
#'
#' @param Data The datapoints to be used. Rows are expected to indicate the points and columes their dimensions
#' @param NumNodes The number of nodes to be used to compute the principal graph
#' @param Method The parametriation associated with a specific graph grammar.
#' @param ReduceMethod An optional string indicating if dimensionality reduction via PCA need to per performed
#' @param ReduceSize Th number of dimensions to retain in case of dimensionality reduction
#' @param ... Additional parameters to be passed to the internal function thqt computes the principal graph
#'
#' @return A list describing the principal elastic graph
#' @export
#'
#' @examples
computeElasticPrincipalGraph <- function(Data, NumNodes, Method = NULL,
                                         ReduceMethod = 'none', ReduceSize = 3, ...) {


  ExtraArgs <- list(...)

  # Reduce dimensions?

  RedDone <- FALSE

  if(ReduceMethod != 'none'){

    print(paste("Performing dimensionality reduction using", ReduceMethod))

    if(ReduceMethod == 'irlba-Lanczos' & ReduceSize == ncol(Data)){

      print(paste("Performing dimensionality reduction using base-svd instead of irlba-Lanczos (RTM)"))

    } else {

      print(paste("Performing dimensionality reduction using", ReduceMethod))

    }

    FunArg <- append(list(Data = Data, Components = ReduceSize, Method = ReduceMethod), ExtraArgs[OkArgs])
    RedData <- do.call(SelectComputePCA, FunArg)

    NewData <- t(t(Data)-RetVal$centers)%*%RetVal$Comp


    if((ReduceMethod != 'base-svd' & ReduceMethod != 'base-cov') & ReduceSize < ncol(Data)){
      print(paste(signif(100*sum(RetVal$ExpVar), 4), "% (estimated) of the original variance vas retained", sep=''))
    } else {
      print(paste(signif(100*sum(RetVal$ExpVar), 4), "% of the original variance vas retained", sep=''))
    }

  } else {

    NewData <- Data
    RedDone <- TRUE

  }


  if(!RedDone){

    print('Invalid dimensionality reduction method selected')
    print('Original data will be used')

  }

  if(is.null(Method)){
    print('Please select the correct parametrization')
  } else {
    Parameters <- do.call(Method, list())
    PrintGraph <- do.call(computeElPT,
                          append(list(Data = NewData, NumNodes = NumNodes, Parameters = Parameters),
                                 ExtraArgs))
    return(PrintGraph)

  }

}

