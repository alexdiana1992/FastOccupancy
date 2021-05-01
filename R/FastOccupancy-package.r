#' Fast fitting of Bayesian occupancy models
#'
#' The package implements the methods to fit Bayesian occupancy models as
#' described in Diana et al. (2021). The model performs estimate of temporal and
#' spatial trends in the occupancy probability and can incorporate covariates.
#' The inference is performed using the  Polya-Gamma scheme for efficient inference.
#' 
#' @details 
#' 
#' Example of the use of the package can be found in the vignette.
#' 
#' @references Diana A., Dennis E., Matechou E. and B.J.T. Morgan. 
#' Fast Bayesian inference for large occupancy datasets, using the Polya-Gamma scheme
#' 
#' @useDynLib FastOccupancy
#' @importFrom Rcpp sourceCpp
#'
#' @examples
#' modelResults <- runModel(sampleData, 
#'                          index_year = 1, 
#'                          index_site = 2, 
#'                          index_occ = 8, 
#'                          index_spatial_x = 3, 
#'                          index_spatial_y = 4, 
#'                          covariates_psi_text = "5", 
#'                          covariates_p_text = "6-7", 
#'                          usingSpatial = TRUE,
#'                          gridStep = .2, 
#'                          nchain = 1, 
#'                          nburn = 100,
#'                          niter = 100) 
#'                          
#' plotOccupancyIndex(modelResults)
#'
#' @docType package
#' @name FastOccupancy
NULL
#> NULL