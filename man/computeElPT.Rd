% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/principal_graphs.R
\name{computeElPT}
\alias{computeElPT}
\title{Compute the the elastic principal graph}
\usage{
computeElPT(Data, NumNodes, Parameters, NodesPositions = NULL, Edges = NULL,
  ...)
}
\arguments{
\item{Data}{The data points to be used to compute the graph}

\item{NumNodes}{The number of vertices to be used to construct the graph}

\item{Parameters}{The name of a parameter function}

\item{...}{Additional parameter to be procesed by the Java procedure. Currently only EP, RP and TrimRadius are supported}
}
\value{
A list describing the principal elastic graph
}
\description{
The function depends on a working JVM with the VDAOEngine.jar library in the libpath.
This is an internal function and should not be accessed directly by the user.
}
