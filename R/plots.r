#' Yearly occupancy probabilities
#' 
#' @description Plots the estimates of the average yearly occupancy
#' probabilities across all sites, as described in Diana et al. (2021).
#' 
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A plot with of the estimates. The interval represents the 
#' 95\% credible interval.
#' 
#' @examples
#' 
#' plotYearlyAveragePresenceProb(sampleResults)
#' 
plotYearlyAveragePresenceProb <- function(modelResults){
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  psi_mean_output <- modelResults$modelOutput$psi_mean_output
  psi_mean_output <- apply(psi_mean_output, 3, c)
  
  CI_psi_mean <- apply(psi_mean_output, 2, function(x){
    c(stats::quantile(x, probs = c(0.025, 0.975)), mean(x))
  })
  
  ggplot2::ggplot(data = NULL, ggplot2::aes(x = years,
                                            y = CI_psi_mean[3,],
                                            ymin = CI_psi_mean[1,],
                                            ymax = CI_psi_mean[2,])) + ggplot2::geom_point() + ggplot2::geom_errorbar() + 
    ggplot2::geom_line() + 
    ggplot2::ylab("Occupancy probability") + ggplot2::scale_x_continuous(name = "Year", breaks = years) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
                   axis.title = ggplot2::element_text(size = 20, face = "bold"),
                   axis.text = ggplot2::element_text(size = 13, face = "bold", angle = 90),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.15),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) + 
    ggplot2::ylim(c(0,1))
  
}

#' Estimates of covariates on occupancy probability
#' 
#' @description Plots the estimates of the covariates on the occupancy probability.
#' 
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A plot with of the estimates. The interval represents the 
#' 95\% credible interval.
#' 
#' @examples
#' 
#' plotOccupancyCovariatesEffect(sampleResults)
#' 
plotOccupancyCovariatesEffect <- function(modelResults) {
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
  
  ggplot2::ggplot(data = NULL,
                  ggplot2::aes(
                    x = colnames(betaeffect),
                    y = betaeffect[2, ],
                    ymin = betaeffect[1, ],
                    ymax = betaeffect[3, ]
                  )) + ggplot2::geom_point() + ggplot2::geom_errorbar() +
    ggplot2::ylab("Effect") + ggplot2::scale_x_discrete(name = "Covariates") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
      axis.title = ggplot2::element_text(size = 20, face = "bold"),
      axis.text = ggplot2::element_text(size = 16, face = "bold"),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black")
    )
  
  
}

#' Estimates of detection probability
#' 
#' @description Plots the estimates of the detection probability. This will plot
#' one value if \code{usingYearDetProb} is set to \code{FALSE} and if set to \code{TRUE}
#' will plot as many values as the number of years.
#' 
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A plot with of the estimates. The interval represents the 
#' 95\% credible interval.
#' 
#' @examples
#' 
#' plotBaseLineCaptureProb(sampleResults)
#' 
plotBaseLineCaptureProb <- function(modelResults) {
  
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
    x_axis <- years
  } else {
    x_axis <- ""
  }
  
  p_plot <- ggplot2::ggplot(data = NULL, ggplot2::aes(
    x = x_axis,
    y = peffect[2,],
    ymin = peffect[1,],
    ymax = peffect[3,]
  )) + ggplot2::geom_point() + ggplot2::geom_errorbar() +
    # ylim(c(0,1)) +
    ggplot2::ylab("Value") + 
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
      axis.title = ggplot2::element_text(size = 20, face = "bold"),
      axis.text = ggplot2::element_text(size = 13, face = "bold", angle = 90),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black")
    ) + ggplot2::ylim(c(min(peffect) / 2, max(peffect) * 1.2))
  
  if(usingYearDetProb){
    p_plot <- p_plot + ggplot2::scale_x_continuous(name = "Years", 
                                                   breaks = x_axis)
  } else {
    p_plot <- p_plot + ggplot2::scale_x_discrete(name = "Detection probability")
  }
  
  p_plot
  
}

#' Estimates of covariates on detection probability
#' 
#' @description Plots the estimates of the covariates on the detection probability.
#' 
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A plot with of the estimates. The interval represents the 
#' 95\% credible interval.
#' 
#' @examples
#' 
#' plotDetectionCovariatesEffect(sampleResults)
#' 
plotDetectionCovariatesEffect <- function(modelResults){
  
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
    
    ggplot2::ggplot(data = NULL,
                    ggplot2::aes(
                      x = namesCovariates,
                      y = betaeffect[2, ],
                      ymin = betaeffect[1, ],
                      ymax = betaeffect[3, ]
                    )) + ggplot2::geom_point() + ggplot2::geom_errorbar() +
      ggplot2::ylab("Effect") + ggplot2::scale_x_discrete(name = "Covariates") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
        axis.title = ggplot2::element_text(size = 20, face = "bold"),
        axis.text = ggplot2::element_text(size = 16, face = "bold"),
        panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
        panel.background = ggplot2::element_rect(fill = "white", color = "black")
      )
  } else {
    print("No detection covariates")
  }
  
}