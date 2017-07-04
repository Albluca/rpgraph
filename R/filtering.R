#' Title
#'
#' @param DataMatrix 
#' @param PCAMat 
#' @param Dimensions 
#' @param PCAMode 
#' @param OutThr 
#'
#' @return
#' @export
#'
#' @examples
PCAOutLiers <- function(DataMatrix, PCAMat = NULL, Dimensions, PCAMode = 'base-svd', OutThr = 5) {
  
  if(PCAMode == 'base-svd' & is.null(PCAMat)){
    
    PCAData <- prcomp(DataMatrix, retx = TRUE, center=FALSE, scale.=FALSE)
    PCAMat <- PCAData$x

  }
  
  Centroid <- colMeans(PCAMat[,1:Dimensions])
  SelDists <- fields::rdist(PCAMat[,1:Dimensions], t(as.matrix(Centroid)))
  return(scater::isOutlier(SelDists, nmads = OutThr))
  
}