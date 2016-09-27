

# PCA Functions --------------------------------------------


#' Compute the PCA for a given numeric matrix usig different methodologies
#'
#' @param DataMatrix The matrix used to compute PCA.
#' @param Components The number of principal components to be retained.
#' @param Method The method to be used to compute PCA. The current implementation accepts 'base-svd',
#' 'base-cov', 'flashPCA-eigen', 'flashPCA-svd', and 'irlba-Lanczos'
#' @return A list composed five elements rotation
#' @return Comp the rotation matrix
#' @return Center A vector containing the velue used to center the original matrix
#'
#'
#' @examples
#' add(1, 1)
#' add(10, 1)

SelectComputePCA <- function(DataMatrix, Components = NULL, Method = 'base-svd', ...){

  tic()

  if(is.null(Components)){
    Components <- min(c(nrow(DataMatrix), ncol(DataMatrix)))/10
    Components <- ceiling(Components)
  }

  ExtraArgs <- list(...)

  if(ncol(DataMatrix) > nrow(DataMatrix)){

    print("The number of columns is largr than the number of rows.")
    print("This is incompatible with some PCA functions")

  }

  if(!is.numeric(Components)){

    print(paste("Invalid number of components, using", ncol(DataMatrix)))
    Components = ncol(DataMatrix)

  } else {

    if((Components > ncol(DataMatrix)) | (Components < 1)){

      print(paste("Invalid number of components, using", ncol(DataMatrix)))
      Components = ncol(DataMatrix)

    }
  }

  if(Method == 'base-svd'){
    # Using svd decomposition
    print("Using base R svd PCA")
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(prcomp, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(x = DataMatrix)
    }

    if(!any(names(FunArg)=="center")){
      FunArg <- append(FunArg, list(center=TRUE))
    }

    PCA <- do.call(what = prcomp, args = FunArg)

    RetVal <- list()
    RetVal$centers <- PCA$center
    RetVal$Comp <- PCA$rotation[,1:Components]
    RetVal$ExpVar <- (PCA$sdev[1:Components]^2)/sum(PCA$sdev^2)

    toc()
    return(RetVal)
  }


  if(Method == 'base-cov'){
    # Using svd decomposition
    print("Using base R covariance PCA")
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(princomp, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(x = DataMatrix)
    }

    if(any(names(FunArg)=="scores")){
      FunArg[["scores"]] <- TRUE
      print("The scores argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(scores=TRUE))
    }

    PCA <- do.call(what = princomp, args = FunArg)

    RetVal <- list()
    RetVal$Comp <- PCA$loadings[,1:Components]
    RetVal$centers <- PCA$center
    RetVal$ExpVar <- (PCA$sdev[1:Components]^2)/sum(PCA$sdev^2)

    toc()
    return(RetVal)
  }


  if(Method == 'flashPCA-eigen'){
    # Using covariance
    print("Using flashpcaR covariance PCA")
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(X = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(X = DataMatrix)
    }

    if(any(names(FunArg)=="method")){
      FunArg[["method"]] <- "eigen"
      print("The method argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(method="eigen"))
    }

    if(any(names(FunArg)=="ndim")){
      FunArg[["ndim"]] <- Components
      print("The ndim argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(ndim=Components))
    }

    if(any(names(FunArg)=="do_loadings")){
      FunArg[["do_loadings"]] <- TRUE
      print("The do_loadings argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(do_loadings=TRUE))
    }



    if(!any(names(FunArg)=="stand")){
      FunArg <- append(FunArg, list(stand="center"))
    }

    if(!any(names(FunArg)=="return_scale")){
      FunArg <- append(FunArg, list(return_scale=TRUE))
    }

    PCA <- do.call(what = flashpca, args = FunArg)

    RetVal <- list()
    RetVal$Comp <- PCA$loadings
    RetVal$centers <- PCA$center

    # Estimating explained variance

    RetVal$ExpVar <- estimate.eig.vpcs(eigenv = PCA$values, M=DataMatrix, estimated = FALSE)$variance.pcs

    toc()
    return(RetVal)
  }


  if(Method == 'flashPCA-svd'){
    # Using svd decomposition
    print("Using flashpcaR svd PCA")
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(X = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(X = DataMatrix)
    }

    if(any(names(FunArg)=="method")){
      FunArg[["method"]] <- "svd"
      print("The method argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(method="svd"))
    }

    if(any(names(FunArg)=="ndim")){
      FunArg[["ndim"]] <- Components
      print("The ndim argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(ndim=Components))
    }

    if(any(names(FunArg)=="do_loadings")){
      FunArg[["do_loadings"]] <- TRUE
      print("The do_loadings argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(do_loadings=TRUE))
    }



    if(!any(names(FunArg)=="stand")){
      FunArg <- append(FunArg, list(stand="center"))
    }

    if(!any(names(FunArg)=="return_scale")){
      FunArg <- append(FunArg, list(return_scale=TRUE))
    }

    PCA <- do.call(what = flashpca, args = FunArg)

    RetVal <- list()
    RetVal$Comp <- PCA$loadings
    RetVal$centers <- PCA$center

    RetVal$ExpVar <- estimate.eig.vpcs(eigenv = PCA$values, M=DataMatrix, estimated = FALSE)$variance.pcs

    toc()
    return(RetVal)
  }


  if(Method == 'irlba-Lanczos'){
    # Using svd decomposition
    print("Using irlba PCA")
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(x = DataMatrix)
    }

    if(any(names(FunArg)=="n")){
      FunArg[["n"]] <- Components
      print("The n argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(n=Components))
    }


    if(!any(names(FunArg)=="center")){
      FunArg <- append(FunArg, list(center=TRUE))
    }

    PCA <- do.call(what = prcomp_irlba, args = FunArg)

    RetVal <- list()
    RetVal$Comp <- PCA$rotation
    RetVal$centers <- PCA$center

    # Estimating explained variance

    RetVal$ExpVar <- estimate.eig.vpcs(eigenv = PCA$sdev^2, M=DataMatrix, estimated = FALSE)$variance.pcs

    toc()
    return(RetVal)
  }


  if(Method == 'nsprcomp'){
    # Using svd decomposition
    print("Using constrained PCA")
    print("Working")

    # get function arguments

    if(length(ExtraArgs)>0){
      OkArgs <- CheckArguments(nsprcomp, names(ExtraArgs), Filter = FALSE)
      FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
    } else {
      FunArg <- list(x = DataMatrix)
    }

    if(any(names(FunArg)=="ncomp")){
      FunArg[["ncomp"]] <- Components
      print("The ncomp argument has been ovrewritten")
    } else {
      FunArg <- append(FunArg, list(ncomp=Components))
    }


    if(!any(names(FunArg)=="center")){
      FunArg <- append(FunArg, list(center=TRUE))
    }

    PCA <- do.call(what = nsprcomp, args = FunArg)

    RetVal <- list()
    RetVal$Comp <- PCA$rotation
    RetVal$centers <- PCA$center

    # Estimating explained variance

    RetVal$ExpVar <- estimate.eig.vpcs(eigenv = PCA$sdev^2, M=DataMatrix, estimated = FALSE)$variance.pcs

    toc()
    return(RetVal)
  }

  print("Unknown methods. Please check.")
  return(NULL)

}
