#' Title
#'
#' @param DataPath 
#' @param GeneInfoPath 
#' @param CellInfoPath 
#' @param GeneCol 
#' @param CellCol 
#'
#' @return
#' @export
#'
#' @examples
LoadMatrixMarket <- function(DataPath, GeneInfoPath=NULL, CellInfoPath=NULL, GeneCol=1, CellCol=1) {
  
  MatrixData <- Matrix::readMM(file = DataPath)
  
  if(!is.null(GeneInfoPath)){
    
    GeneInfo <- readr::read_delim(file = GeneInfoPath, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
    GeneNames <- unlist(GeneInfo[,GeneCol], use.names = FALSE)
      
    rownames(MatrixData) <- GeneNames
  }
  
  if(!is.null(CellInfoPath)){
    CellInfo <- readr::read_delim(file = CellInfoPath, "\t", escape_double = FALSE, trim_ws = TRUE, col_names = FALSE)
    CellNames <- unlist(CellInfo[,CellCol], use.names = FALSE)
    
    colnames(MatrixData) <- CellNames
  }
  
  return(MatrixData)
  
}