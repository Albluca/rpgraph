

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

  tictoc::tic()

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


  if(Method == 'base-cov'){
    # Using svd decomposition
    print("Using base R covariance PCA")
    warning("Columns will be centered and scaled (by SD)")
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
    RetVal$scale <- PCA$scale
    RetVal$ExpVar <- (PCA$sdev[1:Components]^2)/sum(PCA$sdev^2)

    tictoc::toc()
    return(RetVal)
  }


  if(Method == 'flashPCA-eigen'){
    # Using covariance decomposition
    print("Try using flashpcaR covariance PCA")

    CheckOK <- TRUE

    if(Components == min(dim(DataMatrix))){
      warnings("All the components will be computed using flashPCA ... it may be better to simply use base-svd")
    }

    if(!requireNamespace("flashpcaR", quietly = TRUE)){
      warnings("flashpcaR package not available, using base-svd")
      Method <- 'base-svd'
      CheckOK <- FALSE
    }

    if(CheckOK){
      print("Working")

      # get function arguments

      if(length(ExtraArgs)>0){
        OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
        FunArg <- append(list(X = DataMatrix), ExtraArgs[OkArgs])
      } else {
        FunArg <- list(X = DataMatrix)
      }

      # make sure that method is set to "eigen"

      if(any(names(FunArg)=="method")){
        FunArg[["method"]] <- "eigen"
        print("The method argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(method="eigen"))
      }

      # make sure that ndim is set to Components

      if(any(names(FunArg)=="ndim")){
        FunArg[["ndim"]] <- Components
        print("The ndim argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(ndim=Components))
      }

      # make sure that do_loadings is set to TRUE

      if(any(names(FunArg)=="do_loadings")){
        FunArg[["do_loadings"]] <- TRUE
        print("The do_loadings argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(do_loadings=TRUE))
      }

      # Unless it has been specified by the user, the algorithm will center and normalize (by sd) the data

      if(!any(names(FunArg)=="stand")){
        warning("Columns will be centered and scaled (by SD)")
        FunArg <- append(FunArg, list(stand="sd"))
      }

      # make sure that return_scale is set to TRUE

      if(any(names(FunArg)=="return_scale")){
        FunArg[["do_loadings"]] <- TRUE
        print("The return_scale argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(return_scale=TRUE))
      }


      PCA <- do.call(what = flashpcaR::flashpca, args = FunArg)

      RetVal <- list()
      RetVal$Comp <- PCA$loadings
      RetVal$centers <- PCA$center

      # Estimating explained variance

      RetVal <- list()
      RetVal$Comp <- PCA$loadings
      RetVal$centers <- PCA$center
      RetVal$scale <- PCA$scale

      if(!requireNamespace("bigpca", quietly = TRUE)){

        # If bigpca::estimate.eig.vpcs is not available, the expalined variance will be computed using only
        # the eigenvalue computed. This resuls in an overestimations.

        warnings("bigpca package not available, using explained variance will not be estimated.")
        RetVal$ExpVar <- (PCA$value)/sum(PCA$value)

      } else {

        RetVal$ExpVar <- bigpca::estimate.eig.vpcs(eigenv = PCA$values, M=DataMatrix, estimated = FALSE)$variance.pcs

      }

      tictoc::toc()
      return(RetVal)


    }



  }


  if(Method == 'flashPCA-svd'){
    # Using svd decomposition
    print("Try using flashpcaR svd PCA")

    CheckOK <- TRUE

    if(Components == min(dim(DataMatrix))){
      warnings("All the components will be computed using flashPCA ... it may be better to simply use base-svd")
    }

    if(!requireNamespace("flashpcaR", quietly = TRUE)){
      warnings("flashpcaR package not available, using base-svd")
      Method <- 'base-svd'
      CheckOK <- FALSE
    }

    if(CheckOK){
      print("Working")

      # get function arguments and che if everything is ok

      if(length(ExtraArgs)>0){
        OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
        FunArg <- append(list(X = DataMatrix), ExtraArgs[OkArgs])
      } else {
        FunArg <- list(X = DataMatrix)
      }

      # make sure that method is set to "svd"

      if(any(names(FunArg)=="method")){
        FunArg[["method"]] <- "svd"
        print("The method argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(method="svd"))
      }

      # make sure that ndim is set to Components

      if(any(names(FunArg)=="ndim")){
        FunArg[["ndim"]] <- Components
        print("The ndim argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(ndim=Components))
      }

      # make sure that do_loadings is set to TRUE

      if(any(names(FunArg)=="do_loadings")){
        FunArg[["do_loadings"]] <- TRUE
        print("The do_loadings argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(do_loadings=TRUE))
      }


      # Unless it has been specified by the user, the algorithm will center and normalize (by sd) the data

      if(!any(names(FunArg)=="stand")){
        warning("Columns will be centered and scaled (by SD)")
        FunArg <- append(FunArg, list(stand="sd"))
      }

      # make sure that return_scale is set to TRUE

      if(any(names(FunArg)=="return_scale")){
        FunArg[["do_loadings"]] <- TRUE
        print("The return_scale argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(return_scale=TRUE))
      }

      PCA <- do.call(what = flashpcaR::flashpca, args = FunArg)

      RetVal <- list()
      RetVal$Comp <- PCA$loadings
      RetVal$centers <- PCA$center
      RetVal$scale <- PCA$scale

      if(!requireNamespace("bigpca", quietly = TRUE)){

        # If bigpca::estimate.eig.vpcs is not available, the expalined variance will be computed using only
        # the eigenvalue computed. This resuls in an overestimations.

        warnings("bigpca package not available, using explained variance will not be estimated.")
        RetVal$ExpVar <- (PCA$value)/sum(PCA$value)

      } else {

        RetVal$ExpVar <- bigpca::estimate.eig.vpcs(eigenv = PCA$values, M=DataMatrix, estimated = FALSE)$variance.pcs

      }

      tictoc::toc()
      return(RetVal)
    }

  }


  if(Method == 'irlba-Lanczos'){
    # Using svd decomposition
    print("Try using Using irlba PCA")

    CheckOK <- TRUE

    if(Components == min(dim(DataMatrix))){
      warnings("All the components will be computed using base-svd")
      Method <- 'base-svd'
      CheckOK <- FALSE
    }

    if(!requireNamespace("irlba", quietly = TRUE)){
      warnings("irlba package not available, using base-svd")
      Method <- 'base-svd'
      CheckOK <- FALSE
    }

    if(CheckOK){

      print("Working")

      # get function arguments

      if(length(ExtraArgs)>0){
        OkArgs <- CheckArguments(flashpca, names(ExtraArgs), Filter = FALSE)
        FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
      } else {
        FunArg <- list(x = DataMatrix)
      }

      # make sure that n is set to Components

      if(any(names(FunArg)=="n")){
        FunArg[["n"]] <- Components
        print("The n argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(n=Components))
      }


      # Unless it has been specified by the user, the algorithm will center the data

      if(!any(names(FunArg)=="center")){
        warning("Columns will be centered")
        FunArg <- append(FunArg, list(center=TRUE))
      }

      # Unless it has been specified by the user, the algorithm will scale the data

      if(!any(names(FunArg)=="scale.")){
        warning("Columns will be scaled (by SD)")
        FunArg <- append(FunArg, list(scale.=TRUE))
      }

      PCA <- do.call(what = irlba::prcomp_irlba, args = FunArg)

      RetVal <- list()
      RetVal$Comp <- PCA$rotation
      RetVal$centers <- PCA$center
      RetVal$scale <- PCA$scale

      # Estimating explained variance

      if(!requireNamespace("bigpca", quietly = TRUE)){

        # If bigpca::estimate.eig.vpcs is not available, the expalined variance will be computed using only
        # the eigenvalue computed. This resuls in an overestimations.

        warnings("bigpca package not available, using explained variance will not be estimated.")
        RetVal$ExpVar <- (PCA$sdev^2)/sum(PCA$sdev^2)

      } else {

        RetVal$ExpVar <- bigpca::estimate.eig.vpcs(eigenv = PCA$sdev^2, M=DataMatrix, estimated = FALSE)$variance.pcs

      }

      tictoc::toc()
      return(RetVal)

    }


  }


  if(Method == 'nsprcomp'){
    # Using constrained PCA
    print("Try using constrained PCA")

    CheckOK <- TRUE

    if(Components == min(dim(DataMatrix))){
      warnings("All the components will be computed using nsprcomp. This may take a very long time ...")
    }

    if(!requireNamespace("nsprcomp", quietly = TRUE)){
      warnings("nsprcomp package not available, using base-svd")
      Method <- 'base-svd'
      CheckOK <- FALSE
    }

    if(CheckOK){

      print("Working")

      # get function arguments

      if(length(ExtraArgs)>0){
        OkArgs <- CheckArguments(nsprcomp, names(ExtraArgs), Filter = FALSE)
        FunArg <- append(list(x = DataMatrix), ExtraArgs[OkArgs])
      } else {
        FunArg <- list(x = DataMatrix)
      }

      # make sure that ncomp is set to Components

      if(any(names(FunArg)=="ncomp")){
        FunArg[["ncomp"]] <- Components
        print("The ncomp argument has been ovrewritten")
      } else {
        FunArg <- append(FunArg, list(ncomp=Components))
      }


      # Unless it has been specified by the user, the algorithm will center the data

      if(!any(names(FunArg)=="center")){
        warning("Columns will be centered")
        FunArg <- append(FunArg, list(center=TRUE))
      }

      # Unless it has been specified by the user, the algorithm will scale the data

      if(!any(names(FunArg)=="scale.")){
        warning("Columns will be scaled (by SD)")
        FunArg <- append(FunArg, list(scale.=TRUE))
      }

      PCA <- do.call(what = nsprcomp::nsprcomp, args = FunArg)

      RetVal <- list()
      RetVal$Comp <- PCA$rotation
      RetVal$centers <- PCA$center
      RetVal$scale <- PCA$scale

      # Estimating explained variance

      if(!requireNamespace("bigpca", quietly = TRUE)){

        # If bigpca::estimate.eig.vpcs is not available, the expalined variance will be computed using only
        # the eigenvalue computed. This resuls in an overestimations.

        warnings("bigpca package not available, using explained variance will not be estimated.")
        RetVal$ExpVar <- (PCA$sdev^2)/sum(PCA$sdev^2)

      } else {

        RetVal$ExpVar <- bigpca::estimate.eig.vpcs(eigenv = PCA$sdev^2, M=DataMatrix, estimated = FALSE)$variance.pcs

      }

      tictoc::toc()
      return(RetVal)
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

    # Unless it has been specified by the user, the algorithm will center the data

    if(!any(names(FunArg)=="center")){
      warning("Columns will be centered")
      FunArg <- append(FunArg, list(center=TRUE))
    }

    # Unless it has been specified by the user, the algorithm will scale the data

    if(!any(names(FunArg)=="scale.")){
      warning("Columns will be scaled (by SD)")
      FunArg <- append(FunArg, list(scale.=TRUE))
    }


    PCA <- do.call(what = prcomp, args = FunArg)

    RetVal <- list()
    RetVal$centers <- PCA$center
    RetVal$scale <- PCA$scale
    RetVal$Comp <- PCA$rotation[,1:Components]
    RetVal$ExpVar <- (PCA$sdev[1:Components]^2)/sum(PCA$sdev^2)

    tictoc::toc()
    return(RetVal)

  }


  warning("Unknown methods. Please check.")
  return(NULL)

}
