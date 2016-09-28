# install.packages(c("roxygen2", "testthat", "knitr", "rstudioapi"))
# devtools::install_github("hadley/devtools")

# install.packages(c("formatR", "lintr"))

#
# devtools::use_package("rJava", "Depends")
#
# devtools::use_package("igraph", "Imports")
# devtools::use_package("tictoc", "Imports")
# devtools::use_package("rgl", "Imports")
#
# devtools::use_package("bigpca", "Suggests")
# devtools::use_package("flashpcaR", "Suggests")
# devtools::use_package("irlba", "Suggests")
# devtools::use_package("nsprcomp", "Suggests")
# devtools::use_package("repmis", "Suggests")
# devtools::use_package("plotly", "Suggests")
# devtools::use_package("devtools", "Suggests")
# devtools::use_package("repmis", "Suggests")
# devtools::use_package("downloader", "Suggests")
#
# devtools::use_package("pcaMethods", "Suggests")
#
#
# library(downloader)
#
# UrlLoc <- "https://raw.githubusercontent.com/auranic/Elastic-principal-graphs/master/ElasticPrincipalGraphs_matlab/test_data/circle/simple_circle.data"
# TempLoc <- tempfile(pattern = "simple_circle", fileext = ".data")
# download.file(url = UrlLoc, destfile = TempLoc, mode = "w")
# simple_circle <- read.delim(TempLoc, header=FALSE, stringsAsFactors=FALSE)
# file.remove(TempLoc)
#
# devtools::use_data(simple_circle)
#
#
#
# UrlLoc <- "https://raw.githubusercontent.com/auranic/Elastic-principal-graphs/master/ElasticPrincipalGraphs_matlab/test_data/tree23/tree23.data"
# TempLoc <- tempfile(pattern = "tree23", fileext = ".data")
# download.file(url = UrlLoc, destfile = TempLoc, mode = "w")
# simple_tree <- read.delim(TempLoc, header=FALSE, stringsAsFactors=FALSE)
# file.remove(TempLoc)
#
# devtools::use_data(simple_tree)
#
#
#
# UrlLoc <- "https://raw.githubusercontent.com/auranic/Elastic-principal-graphs/master/ElasticPrincipalGraphs_matlab/test_data/iris/iris.txt"
# TempLoc <- tempfile(pattern = "iris", fileext = ".txt")
# download.file(url = UrlLoc, destfile = TempLoc, mode = "w")
# simple_iris <- read.delim(TempLoc, header=TRUE, stringsAsFactors=FALSE)
# file.remove(TempLoc)
#
# devtools::use_data(simple_iris)

# devtools::use_vignette("Simple_Examples")
