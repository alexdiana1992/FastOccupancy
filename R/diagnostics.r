
#' Traceplot of the occupancy index rate
#' 
#' @description Traceplot of the index occupancy rate, with 
#' effective sample size.
#' 
#' @param modelResults Output of the function \code{runModel}
#' @param index_year Index of the year to plot. Indexes go from 1 to the number of years.
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return The traceplot with the estimate of the effective sample size.
#' 
#' @examples
#' 
#' tracePlot_OccupancyIndexRate(sampleResults, 1)
#' 
tracePlot_OccupancyIndexRate <- function(modelResults, index_year) {
  
  psi_mean_output <- modelResults$modelOutput$psi_mean_output
  
  nchain <- nrow(psi_mean_output)
  niter <- ncol(psi_mean_output)
  
  psi_mean_output <-
    matrix(psi_mean_output[, , index_year], nrow = nchain, ncol = niter)
  
  psi_mean_output_long <- reshape2::melt(psi_mean_output)
  
  diagnosticsPlot <-
    ggplot2::ggplot(data = psi_mean_output_long, ggplot2::aes(
      x = Var2,
      y = value,
      group = Var1,
      color = factor(Var1)
    )) + ggplot2::geom_line() +
    ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 1
      ),
      axis.line = ggplot2::element_line(colour = "black", size = 0.15),
      # panel.grid.minor = element_line(colour="grey", size=0.15),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      legend.position = "none"
    )
  
  
  plotTitle <- createPlotTitle(psi_mean_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
  
  diagnosticsPlot
}

#' Traceplot of the year-specific random effects
#' 
#' @description Traceplot of the year-specific random effects, with 
#' effective sample size.
#' 
#' @param modelResults Output of the function \code{runModel}
#' @param index_year Index of the year of the random effect to plot. Indexes go from 1 to the number of years.
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return The traceplot with the estimate of the effective sample size.
#' 
#' @examples
#' 
#' tracePlot_OccupancyYearEffect(sampleResults, 1)
#' 
tracePlot_OccupancyYearEffect <- function(modelResults, index_year) {
  beta_psi_output <- modelResults$modelOutput$beta_psi_output
  
  nchain <- nrow(beta_psi_output)
  niter <- ncol(beta_psi_output)
  
  beta_psi_output <-
    matrix(beta_psi_output[, , index_year], nrow = nchain, ncol = niter)
  
  beta_psi_output_long <- reshape2::melt(beta_psi_output)
  
  diagnosticsPlot <-
    ggplot2::ggplot(data = beta_psi_output_long, ggplot2::aes(
      x = Var2,
      y = value,
      group = Var1,
      color = factor(Var1)
    )) + ggplot2::geom_line() +
    ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 1
      ),
      axis.line = ggplot2::element_line(colour = "black", size = 0.15),
      # panel.grid.minor = element_line(colour="grey", size=0.15),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      legend.position = "none"
    )
  
  
  plotTitle <- createPlotTitle(beta_psi_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
  
  diagnosticsPlot
}

#' Traceplot of the site-specific autocorrelated random effects
#' 
#' @description Traceplot of the site-specific random effects, with 
#' effective sample size.
#' 
#' @param modelResults Output of the function \code{runModel}
#' @param index_site Index of the random effect to plot. Indexes go from 1 to the number of sites in the approximation.
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return The traceplot with the estimate of the effective sample size.
#' 
#' @examples
#' 
#' tracePlot_OccupancySiteEffect(sampleResults, 1)
#' 
tracePlot_OccupancySiteEffect <- function(modelResults, index_site) {
  usingSpatial <- modelResults$dataCharacteristics$usingSpatial
  
  if (!usingSpatial) {
    print("No spatial effect")
  } else {
    years <- modelResults$dataCharacteristics$Years
    Y <- length(years)
    
    beta_psi_output <- modelResults$modelOutput$beta_psi_output
    
    nchain <- nrow(beta_psi_output)
    niter <- ncol(beta_psi_output)
    
    beta_psi_output <-
      matrix(beta_psi_output[, , Y + index_site], nrow = nchain, ncol = niter)
    
    beta_psi_output_long <- reshape2::melt(beta_psi_output)
    
    diagnosticsPlot <-
      ggplot2::ggplot(data = beta_psi_output_long, ggplot2::aes(
        x = Var2,
        y = value,
        group = Var1,
        color = factor(Var1)
      )) + ggplot2::geom_line() +
      ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
        axis.text.x = ggplot2::element_text(
          size = 11,
          face = "bold",
          hjust = 1
        ),
        axis.line = ggplot2::element_line(colour = "black", size = 0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
        panel.background = ggplot2::element_rect(fill = "white", color = "black"),
        legend.position = "none"
      )
    
    
    plotTitle <- createPlotTitle(beta_psi_output, nchain)
    
    diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
    
    return(diagnosticsPlot)
  }
  
}

#' Traceplot of the covariates of the occupancy probability
#' 
#' @description Traceplot of the covariates coefficient, with 
#' effective sample size.
#' 
#' @param modelResults Output of the function \code{runModel}
#' @param index_cov Index of the covariate to plot. Indexes go from 1 to the number of covariate.
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return The traceplot with the estimate of the effective sample size.
#' 
#' @examples
#' 
#' tracePlot_OccupancyCovariate(sampleResults, 1)
#' 
tracePlot_OccupancyCovariate <- function(modelResults, index_cov) {
  beta_psi_output <- modelResults$modelOutput$beta_psi_output
  
  usingSpatial <- modelResults$dataCharacteristics$usingSpatial
  
  Y <- length(modelResults$dataCharacteristics$Years)
  
  if (usingSpatial) {
    X_centers <- nrow(modelResults$dataCharacteristics$X_tilde)
  } else {
    X_centers <- 0
  }
  
  numTimeSpaceCov <-
    modelResults$dataCharacteristics$numTimeSpaceCov
  
  if ((Y + X_centers + numTimeSpaceCov + index_cov) <= dim(beta_psi_output)[3]) {
    nchain <- nrow(beta_psi_output)
    niter <- ncol(beta_psi_output)
    
    beta_psi_output <-
      matrix(beta_psi_output[, , Y + X_centers + numTimeSpaceCov + index_cov], nrow = nchain, ncol = niter)
    
    beta_psi_output_long <- reshape2::melt(beta_psi_output)
    
    diagnosticsPlot <-
      ggplot2::ggplot(data = beta_psi_output_long, ggplot2::aes(
        x = Var2,
        y = value,
        group = Var1,
        color = factor(Var1)
      )) + ggplot2::geom_line() +
      ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
        axis.text.x = ggplot2::element_text(
          size = 11,
          face = "bold",
          hjust = 1
        ),
        axis.line = ggplot2::element_line(colour = "black", size = 0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
        panel.background = ggplot2::element_rect(fill = "white", color = "black"),
        legend.position = "none"
      )
    
    plotTitle <- createPlotTitle(beta_psi_output, nchain)
    
    diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
    
    return(diagnosticsPlot)
    
  } else {
    print("Index above the number of covariates")
    
  }
  
}

#' Traceplot of the intercept of the detection probability
#' 
#' @description Traceplot of the intercept, with 
#' effective sample size.
#' 
#' @param modelResults Output of the function \code{runModel}
#' @param index_year Index of the year to plot, in case \code{usingYearDetProb} is 
#' set to \code{TRUE}. Indexes go from 1 to the number of years.
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return The traceplot with the estimate of the effective sample size.
#' 
#' @examples
#' 
#' tracePlot_DetectionIntercept(sampleResults, 1)
#' 
tracePlot_DetectionIntercept <- function(modelResults, index_year) {
  
  usingYearDetProb <- modelResults$dataCharacteristics$usingYearDetProb
  
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  if(!usingYearDetProb){
    index_year <- 1
  } else if(index_year > Y){
    
    print("Index above the number of years")
    return(NULL)  
    
  }
  
  beta_p_output <- modelResults$modelOutput$beta_p_output
  
  nchain <- nrow(beta_p_output)
  niter <- ncol(beta_p_output)
  
  beta_p_output <-
    matrix(beta_p_output[, , index_year], nrow = nchain, ncol = niter)
  
  beta_psi_output_long <- reshape2::melt(beta_p_output)
  
  diagnosticsPlot <-
    ggplot2::ggplot(data = beta_psi_output_long, ggplot2::aes(
      x = Var2,
      y = value,
      group = Var1,
      color = factor(Var1)
    )) + ggplot2::geom_line() +
    ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
      axis.title = ggplot2::element_text(size = 16, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      axis.text.x = ggplot2::element_text(
        size = 11,
        face = "bold",
        hjust = 1
      ),
      axis.line = ggplot2::element_line(colour = "black", size = 0.15),
      # panel.grid.minor = element_line(colour="grey", size=0.15),
      panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black"),
      legend.position = "none"
    )
  
  plotTitle <- createPlotTitle(beta_p_output, nchain)
  
  diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
  
  diagnosticsPlot
}

#' Traceplot of the covariates of the detection probability
#' 
#' @description Traceplot of the covariates coefficient, with 
#' effective sample size.
#' 
#' @param modelResults Output of the function \code{runModel}
#' @param index_cov Index of the covariate to plot. Indexes go from 1 to the number of covariate.
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return The traceplot with the estimate of the effective sample size.
#' 
#' @examples
#' 
#' tracePlot_DetectionCovariate(sampleResults, 1)
#' 
tracePlot_DetectionCovariate <- function(modelResults, index_cov) {
  beta_p_output <- modelResults$modelOutput$beta_p_output
  
  usingYearDetProb <- modelResults$dataCharacteristics$usingYearDetProb
  p_intercepts <- ifelse(usingYearDetProb, Y, 1)
  
  if ((p_intercepts + index_cov) <= dim(beta_p_output)[3]) {
    nchain <- nrow(beta_p_output)
    niter <- ncol(beta_p_output)
    
    beta_p_output <-
      matrix(beta_p_output[, , p_intercepts + index_cov], nrow = nchain, ncol = niter)
    
    beta_psi_output_long <- reshape2::melt(beta_p_output)
    
    diagnosticsPlot <-
      ggplot2::ggplot(data = beta_psi_output_long, ggplot2::aes(
        x = Var2,
        y = value,
        group = Var1,
        color = factor(Var1)
      )) + ggplot2::geom_line() +
      ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
        axis.title = ggplot2::element_text(size = 16, face = "bold"),
        axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
        axis.text.x = ggplot2::element_text(
          size = 11,
          face = "bold",
          hjust = 1
        ),
        axis.line = ggplot2::element_line(colour = "black", size = 0.15),
        # panel.grid.minor = element_line(colour="grey", size=0.15),
        panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
        panel.background = ggplot2::element_rect(fill = "white", color = "black"),
        legend.position = "none"
      )
    
    plotTitle <- createPlotTitle(beta_p_output, nchain)
    
    diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
    
    return(diagnosticsPlot)
    
  } else {
    print("Index above the number of covariates")
    
  }
  
}

# tracePlot_l_T <- function(modelResults) {
#   l_T_output <- modelResults$modelOutput$l_T_output
#   
#   nchain <- nrow(l_T_output)
#   
#   l_T_psi_output_long <- reshape2::melt(l_T_output)
#   
#   diagnosticsPlot <-
#     ggplot(data = l_T_psi_output_long, aes(
#       x = Var2,
#       y = value,
#       group = Var1,
#       color = factor(Var1)
#     )) + geom_line() +
#     xlab("Iterations") + ylab("Value") +
#     theme(
#       plot.title = element_text(hjust = 0.5, size = 17),
#       axis.title = element_text(size = 16, face = "bold"),
#       axis.text.y = element_text(size = 11, face = "bold"),
#       axis.text.x = element_text(
#         size = 11,
#         face = "bold",
#         hjust = 1
#       ),
#       axis.line = element_line(colour = "black", size = 0.15),
#       # panel.grid.minor = element_line(colour="grey", size=0.15),
#       panel.grid.major = element_line(colour = "grey", size = 0.15),
#       panel.background = element_rect(fill = "white", color = "black"),
#       legend.position = "none"
#     )
#   
#   plotTitle <- createPlotTitle(l_T_output, nchain)
#   
#   diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
#   
#   diagnosticsPlot
#   
# }
# 
# tracePlot_sigma_T <- function(modelResults) {
#   sigma_T_output <- modelResults$modelOutput$sigma_T_output
#   
#   nchain <- nrow(sigma_T_output)
#   
#   sigma_T_output_long <- reshape2::melt(sigma_T_output)
#   
#   diagnosticsPlot <-
#     ggplot(data = sigma_T_output_long, aes(
#       x = Var2,
#       y = value,
#       group = Var1,
#       color = factor(Var1)
#     )) + geom_line() +
#     xlab("Iterations") + ylab("Value") +
#     theme(
#       plot.title = element_text(hjust = 0.5, size = 17),
#       axis.title = element_text(size = 16, face = "bold"),
#       axis.text.y = element_text(size = 11, face = "bold"),
#       axis.text.x = element_text(
#         size = 11,
#         face = "bold",
#         hjust = 1
#       ),
#       axis.line = element_line(colour = "black", size = 0.15),
#       # panel.grid.minor = element_line(colour="grey", size=0.15),
#       panel.grid.major = element_line(colour = "grey", size = 0.15),
#       panel.background = element_rect(fill = "white", color = "black"),
#       legend.position = "none"
#     )
#   
#   plotTitle <- createPlotTitle(sigma_T_output, nchain)
#   
#   diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
#   
#   diagnosticsPlot
#   
# }
# 
# tracePlot_l_S <- function(modelResults) {
#   usingSpatial <- modelResults$dataCharacteristics$usingSpatial
#   
#   if (usingSpatial) {
#     l_s_output <- modelResults$modelOutput$l_s_output
#     
#     nchain <- nrow(l_s_output)
#     
#     l_s_psi_output_long <- reshape2::melt(l_s_output)
#     
#     diagnosticsPlot <-
#       ggplot(data = l_s_psi_output_long, aes(
#         x = Var2,
#         y = value,
#         group = Var1,
#         color = factor(Var1)
#       )) + geom_line() +
#       xlab("Iterations") + ylab("Value") +
#       theme(
#         plot.title = element_text(hjust = 0.5, size = 17),
#         axis.title = element_text(size = 16, face = "bold"),
#         axis.text.y = element_text(size = 11, face = "bold"),
#         axis.text.x = element_text(
#           size = 11,
#           face = "bold",
#           hjust = 1
#         ),
#         axis.line = element_line(colour = "black", size = 0.15),
#         # panel.grid.minor = element_line(colour="grey", size=0.15),
#         panel.grid.major = element_line(colour = "grey", size = 0.15),
#         panel.background = element_rect(fill = "white", color = "black"),
#         legend.position = "none"
#       )
#     
#     plotTitle <- createPlotTitle(l_s_output, nchain)
#     
#     diagnosticsPlot <- diagnosticsPlot + ggtitle(plotTitle)
#     
#     return(diagnosticsPlot)
#   } else {
#     print("No spatial effect")
#   }
#   
#   
# }
# 
# tracePlot_sigma_S <- function(modelResults) {
#   sigma_s_output <- modelResults$modelOutput$sigma_s_output
#   
#   nchain <- nrow(sigma_s_output)
#   
#   sigma_s_output_long <- reshape2::melt(sigma_s_output)
#   
#   diagnosticsPlot <-
#     ggplot2::ggplot(data = sigma_s_output_long, aes(
#       x = Var2,
#       y = value,
#       group = Var1,
#       color = factor(Var1)
#     )) + ggplot2::geom_line() +
#     ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
#     ggplot2::theme(
#       plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
#       axis.title = ggplot2::element_text(size = 16, face = "bold"),
#       axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
#       axis.text.x = ggplot2::element_text(
#         size = 11,
#         face = "bold",
#         hjust = 1
#       ),
#       axis.line = ggplot2::element_line(colour = "black", size = 0.15),
#       # panel.grid.minor = element_line(colour="grey", size=0.15),
#       panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
#       panel.background = ggplot2::element_rect(fill = "white", color = "black"),
#       legend.position = "none"
#     )
#   
#   plotTitle <- createPlotTitle(sigma_s_output, nchain)
#   
#   diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
#   
#   diagnosticsPlot
#   
# }
# 
# tracePlot_sigma_eps <- function(modelResults) {
#   sigma_as_output <- modelResults$modelOutput$sigma_as_output
#   
#   nchain <- nrow(sigma_as_output)
#   
#   sigma_as_output_long <- reshape2::melt(sigma_as_output)
#   
#   ggplot2::diagnosticsPlot <-
#     ggplot2::ggplot(data = sigma_as_output_long, aes(
#       x = Var2,
#       y = value,
#       group = Var1,
#       color = factor(Var1)
#     )) + ggplot2::geom_line() +
#     ggplot2::xlab("Iterations") + ggplot2::ylab("Value") +
#     ggplot2::theme(
#       plot.title = ggplot2::element_text(hjust = 0.5, size = 17),
#       axis.title = ggplot2::element_text(size = 16, face = "bold"),
#       axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
#       axis.text.x = ggplot2::element_text(
#         size = 11,
#         face = "bold",
#         hjust = 1
#       ),
#       axis.line = ggplot2::element_line(colour = "black", size = 0.15),
#       # panel.grid.minor = element_line(colour="grey", size=0.15),
#       panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
#       panel.background = ggplot2::element_rect(fill = "white", color = "black"),
#       legend.position = "none"
#     )
#   
#   plotTitle <- createPlotTitle(sigma_as_output, nchain)
#   
#   diagnosticsPlot <- diagnosticsPlot + ggplot2::ggtitle(plotTitle)
#   
#   diagnosticsPlot
#   
# }

createPlotTitle <- function(mcmc_output, nchain) {
  eff_samplesize <- ess(mcmc_output)
  
  plotTitle <-
    paste0("Effective sample size = ", round(eff_samplesize, 1))
  
  if (nchain > 1) {
    Rhat <- compute_rhat(mcmc_output)
    
    plotTitle <-
      paste0(plotTitle, paste0(" / R.hat = ", round(Rhat, 3)))
  }
  
  plotTitle
}

compute_rhat <- function(mcmc_output) {
  mcmc_output_list <- lapply(1:nrow(mcmc_output), function(i) {
    coda::mcmc(mcmc_output[i, ])
  })
  mcmc_output_list_2 <- coda::as.mcmc.list(mcmc_output_list)
  
  Rhat <- coda::gelman.diag(mcmc_output_list_2)
  
  Rhat$psrf[2]
}

ess <- function(mcmc_output) {
  mcmc_output_list <- lapply(1:nrow(mcmc_output), function(i) {
    coda::mcmc(mcmc_output[i, ])
  })
  mcmc_output_list_2 <- coda::as.mcmc.list(mcmc_output_list)
  
  coda::effectiveSize(mcmc_output_list_2)
  
}