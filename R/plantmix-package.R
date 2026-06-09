#' @import ggplot2
#' @importFrom graphics axis image pairs par title
#' @importFrom grDevices grey.colors
#' @importFrom igraph add_edges as_data_frame degree delete_edges graph_from_adjacency_matrix is_bipartite laplacian_matrix layout_as_bipartite sample_bipartite sample_gnm set_vertex_attr vertex_attr which_multiple
#' @import lme4
#' @importClassesFrom Matrix dgTMatrix
#' @importFrom MASS mvrnorm
#' @importFrom Matrix diag Diagonal drop0 Matrix rowSums sparseMatrix summary t
#' @importFrom methods as is
#' @importFrom Rcpp evalCpp
#' @importFrom stats as.formula coef contr.sum cor extractAIC fitted formula lm model.matrix pnorm qbeta rnorm nlminb residuals setNames
#' @import TMB
#' @importFrom utils capture.output head str strcapture
#' @rawNamespace useDynLib(plantmix, .registration=TRUE); useDynLib(plantmix_TMBExports)
NULL

## see https://github.com/cran/igraph/blob/master/R/igraph-package.R#L28:L36

#' The plantmix package
#'
#' plantmix is a package to study plant mixtures.
#'
#' @rdname aaa-plantmix-package
#' @name plantmix-package
#' @keywords internal
#' @aliases plantmix-package plantmix
#'
#' @section Introduction:
#' The main goal of the plantmix package is to provide a set of functions
#' for pain-free analysis of both varietal (intraspecific) and
#' crop (interspecific) mixtures, such as
#' optimizing experimental designs, simulating data sets, fitting
#' GMA-SMA and DBV-SBV models, computing diversity indices, etc.
#' Work is ongoing on DBV-SBV models for varietal mixtures, spatial heterogeneity, etc.
#'
#' @section Optimizing experimental designs:
#' \itemize{
#'   \item \code{\link{getDesignBinaryVarMix}}, \code{\link{getDesignBinaryCropMix}}
#'   \item \code{\link{plotDesignVarMix}}, \code{\link{plotDesignCropMix}}
#'   \item \code{\link{plotDiallel}}, \code{\link{getMixtureList}}, \code{\link{getMixturesPerGeno}}
#' }
#'
#' @section Fitting GMA-SMA models:
#' \itemize{
#'   \item \code{\link{mkZGMA}}, \code{\link{mkZSMA}}, \code{\link{mkAllZSMA}}
#'   \item \code{\link{fitGMASMA}}, \code{\link{summarizeGMASMAs}}
#' }
#'
#' @section Fitting DBV-SBV models:
#' \itemize{
#'   \item \code{\link{pivotMixData2Long}}, \code{\link{pivotMixData2Wide}}
#'   \item \code{\link{mkZinterspe}}, \code{\link{fitDBVSBVinter}}
#' }
#'
#' @section Simulating mixture data:
#' \itemize{
#'   \item \code{\link{simulDBVSBVinter}}
#' }
#'
#' @section Interaction indices:
#' \itemize{
#'   \item \code{\link{RYM}}
#'   \item \code{\link{RY}}, \code{\link{RYT}}, \code{\link{LER}}, \code{\link{RYP}}
#'   \item \code{\link{CC}}, \code{\link{RII}}, \code{\link{RIInet}}
#'   \item \code{\link{estimRYMRep}}, \code{\link{estimRYRep}}
#' }
#'
#' @section Miscellaneous:
#' \itemize{
#'   \item \code{\link{rmatnorm}}, \code{\link{vec}}, \code{\link{invvec}}, \code{\link{imageMat}}
#'   \item \code{\link{simulGenosDoseStruct}}, \code{\link{estimGRM}}
#'   \item \code{\link{paramBoot4TMB}}
#'   \item \code{\link{mixSowingWeight}}, \code{\link{nbSeedsToSownInPure}}
#' }
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL
