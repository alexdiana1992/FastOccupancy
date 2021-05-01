#' Occupancy index
#' 
#' @description Prints the estimates of the occupancy index, 
#' as described in Diana et al. (2021).
#' 
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A data.frame with of the 95\% credible interval and the median.
#' 
#' @examples
#' 
#' printOccupancyIndex(sampleResults)
#' 
printOccupancyIndex <- function(modelResults){
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  psi_mean_output <- modelResults$modelOutput$psi_mean_output
  psi_mean_output <- apply(psi_mean_output, 3, c)
  
  CI_psi_mean <- apply(psi_mean_output, 2, function(x){
    c(stats::quantile(x, probs = c(0.025, .5, 0.975)))
  })
  
  colnames(CI_psi_mean) <- years
  
  CI_psi_mean <- as.data.frame(CI_psi_mean)
  
  CI_psi_mean
}

#' Estimates of covariates on occupancy probability
#' 
#' @description Prints the estimates of the covariates on the occupancy probability.
#' 
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A data.frame with of the 95\% credible interval and the median.
#' 
#' @examples
#' 
#' printOccupancyCovariatesEffect(sampleResults)
#' 
printOccupancyCovariatesEffect <- function(modelResults) {
  usingSpatial <- modelResults$dataCharacteristics$usingSpatial
  
  X_tilde <- modelResults$dataCharacteristics$X_tilde
  if (usingSpatial) {
    X_centers <- nrow(X_tilde)
  } else {
    X_centers <- 0
  }
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  numTimeSpaceCov <-
    modelResults$dataCharacteristics$numTimeSpaceCov
  
  namesCovariates <-
    modelResults$dataCharacteristics$nameVariables_psi
  ncov_psi <- length(namesCovariates)
  
  if (usingSpatial) {
    namesCovariates_xt <- c("X_T", "Y_T")
  } else {
    namesCovariates_xt <- c()
  }
  
  namesCovariates_all <- c(namesCovariates_xt, namesCovariates)
  ncov_psi_all <- ncov_psi + numTimeSpaceCov
  
  beta_psi_output <- modelResults$modelOutput$beta_psi_output
  
  beta_psi_output <- apply(beta_psi_output, 3, c)
  
  betaeffect <-
    apply(beta_psi_output[, Y + X_centers + 1:ncov_psi_all], 2, function(x) {
      stats::quantile(x, probs = c(0.025, 0.5, 0.975))
    })
  
  colnames(betaeffect) <- namesCovariates_all
  
  betaeffect <- as.data.frame(betaeffect)
  
  betaeffect  
}

#' Estimates of detection probability
#' 
#' @description Prints the estimates of the detection probability. 
#' 
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A data.frame with the 95\% credible interval and the median.
#' 
#' @examples
#' 
#' printBaseLineCaptureProb(sampleResults)
#' 
printBaseLineCaptureProb <- function(modelResults) {
  
  usingYearDetProb <- modelResults$dataCharacteristics$usingYearDetProb
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  p_intercepts <- ifelse(usingYearDetProb, Y, 1)
  
  beta_p_output <- modelResults$modelOutput$beta_p_output
  
  beta_p_output <- apply(beta_p_output, 3, c)
  
  peffect <-
    apply(beta_p_output[, 1:p_intercepts,drop = F], 2, function(x) {
      stats::quantile(logit(x), probs = c(0.025, 0.5, 0.975))
    })
  
  if(usingYearDetProb){
    colnames(peffect) <- years
  } else {
    colnames(peffect) <- ""
  }
  
  peffect <- as.data.frame(peffect)
  
  peffect
}

#' Estimates of covariates on detection probability
#' 
#' @description Prints the estimates of the covariates on the detection probability.
#' 
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A data.frame with the 95\% credible interval and the median.
#' 
#' @examples
#' 
#' printDetectionCovariatesEffect(sampleResults)
#' 
printDetectionCovariatesEffect <- function(modelResults){
  
  usingYearDetProb <- modelResults$dataCharacteristics$usingYearDetProb
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  p_intercepts <- ifelse(usingYearDetProb, Y, 1)
  
  namesCovariates <- modelResults$dataCharacteristics$nameVariables_p
  ncov_p <- length(namesCovariates)
  
  if (ncov_p > 0) {
    beta_p_output <- modelResults$modelOutput$beta_p_output
    
    beta_p_output <- apply(beta_p_output, 3, c)
    
    betaeffect <- apply(beta_p_output[, -(1:p_intercepts), drop = F], 2, function(x) {
      stats::quantile(x, probs = c(0.025, 0.5, 0.975))
    })
    
    betaeffect <- as.data.frame(betaeffect)
    
    return(betaeffect)
  } else {
    print("No detection covariates")
    
    return(NULL)
  }
  
}
