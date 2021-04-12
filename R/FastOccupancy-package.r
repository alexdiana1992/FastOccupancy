#' foo: A package for computating the notorious bar statistic
#'
#' The foo package provides three categories of important functions:
#' foo, bar and baz.
#' 
#' @useDynLib FastOccupancy
#' @importFrom Rcpp sourceCpp
#' 
#' @section FastOccupancy functions:
#' The foo functions ...
#'
#' @examples
#' simulateData(Y = 3, S_years = rep(50, 3),
#'              V_lambda = 2, mu_psi = 0, b_t = c(1,2,3), beta_psi_tsp = c(0,0),
#'              beta_psi = c(1), X = matrix(runif(50 * 2), 50, 2), a_s_site = rnorm(50),
#'               mu_p = -1, beta_p = 0)
#'
#' @docType package
#' @name FastOccupancy
NULL
#> NULL