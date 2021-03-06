#' Count pubmed results
#'
#' @param Tokens 
#'
#' @return
#' @export
#'
#' @examples
CountPubmedResults <- function(Tokens) {
  
  # Tokens <- list("cell cycle", "fos")
  
  SrcStr <- ""
  for(i in 1:length(Tokens)){
    SrcStr <- paste(SrcStr, "+%22", sub(Tokens[[i]], pattern = " ", replacement = "+", fixed = TRUE), '%22', sep = "")
  }
  
  url <- paste("http://www.ncbi.nlm.nih.gov/pubmed/?term=", SrcStr, sep = "")
  
  hh <-  xml2::read_html(url, options = c("RECOVER", "NOERROR", "NOBLANKS"))
  src <- XML::htmlTreeParse(hh,useInternalNodes=TRUE)
  tags <- XML::xpathApply(src, "//meta[@name='ncbi_resultcount']", XML::xmlAttrs)
  return(as.numeric(tags[[1]]["content"]))
  
}