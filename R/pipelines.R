# A function to standardize the Initial steps of the analysis

StudyCellCycle <- function(ExpressionMatrix, GeneDetectedFilter = 2.5, GeneCountFilter = 2.5,
                           MinCellExp = 1, VarFilter = 0, LogTranform = TRUE, Centering = FALSE, Scaling = FALSE) {
  
  print(paste("Expression matrix contains", nrow(ExpressionMatrix), "cells and", nrow(ExpressionMatrix), "genes"))
  
  # Filtering ---------------------------------------------------------------
  
  print("Stage I - Cell filtering")
  
  hist(apply(ExpressionMatrix > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
       freq = TRUE, ylab = "Number of cells")
  
  OutExpr <- scater::isOutlier(rowSums(ExpressionMatrix>0), nmads = GeneDetectedFilter)
  
 
  
  
  hist(apply(ExpressionMatrix, 1, sum), main = "Reads per cell", xlab = "Reads count",
       freq = TRUE, ylab = "Number of cells")
  
  OutCount <- scater::isOutlier(rowSums(ExpressionMatrix>0), nmads = GeneCountFilter)

  
  print(paste(sum(OutExpr & OutCount), "Cells will be removed due to both filtering"))
  print(paste(sum(OutExpr), "Cells will be removed due to gene count filtering"))
  print(paste(sum(OutCount), "Cells will be removed due to read count filtering"))
  print(paste(sum(!(OutExpr & OutCount)), "Cells will be used for analysis")) 

  NormExpressionMatrix <- ExpressionMatrix[!(OutExpr & OutCount),]
  
  readline("Press any key")
  
  par(mfcol=c(1,2))
  
  hist(apply(NormExpressionMatrix > 0, 1, sum), main = "Genes per cell", xlab = "Genes count",
       freq = TRUE, ylab = "Number of cells (After filtering)")
  
  hist(apply(NormExpressionMatrix, 1, sum), main = "Reads per cell", xlab = "Reads count",
       freq = TRUE, ylab = "Number of cells (After filtering)")
  
  
  print("Stage II - Genes filtering")
  
  print(paste("Removing genes that are detected in less than", MinCellExp, "cell(s)"))
  
  GeneFilterBool <- apply(NormExpressionMatrix > 0, 2, sum) < MinCellExp
  NormExpressionMatrix <- NormExpressionMatrix[, !GeneFilterBool]

  print(paste(sum(GeneFilterBool), "genes will be removed"))
  print(paste(ncol(NormExpressionMatrix), "genes will be used for the analysis"))

  if(LogTranform){
    print("Transforming using pseudo counts (Log10(x+1))")
    NormExpressionMatrix <- log10(NormExpressionMatrix + 1)
  }
  
  print("Plotting variance distribution (For reference only ATM)")
    
  par(mfcol=c(1,2))

  VarVect <- apply(NormExpressionMatrix, 2, var)
  
  plot(density(VarVect, from=min(VarVect)), xlab = "Variance", main = "Variance distribution")
  abline(v=quantile(VarVect, c(.01, .05, .25, .50)), lty=2, col=c('red', 'green', 'blue', "black"))
  legend(x = "topright", legend = c("1% Q", "5% Q", "25% Q", "50% Q"),
         col=c('red', 'green', 'blue', "black"), lty=2)
  
  plot(density(VarVect, from=min(VarVect)), xlab = "Variance (1% - 75% Q)",
       xlim = quantile(VarVect, c(.01, .75)), main = "Variance distribution")
  abline(v=quantile(VarVect, c(.01, .05, .25, .50)), lty=2, col=c('red', 'green', 'blue', "black"))
    
  legend(x = "topright", legend = c("1% Q", "5% Q", "25% Q", "50% Q"),
         col=c('red', 'green', 'blue', "black"), lty=2)
  
  print("Stage II - Gene scaling")
  
  if(Centering){
    print("Gene expression will be centered")
  } 
  
  if(Scaling){
    print("Gene expression will be scaled")
  } 
  
  NormExpressionMatrix <- scale(NormExpressionMatrix, Centering, Scaling)
  
  print("Stage II - PCA")
  
  
}


