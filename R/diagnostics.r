
#' Traceplot of the year-specific random effects
#' 
#' @description Traceplot of the year-specific random effects, with 
#' effective sample size..
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

tracePlot_DetectionIntercept <- function(modelResults) {
  beta_p_output <- modelResults$modelOutput$beta_p_output
  
  nchain <- nrow(beta_p_output)
  niter <- ncol(beta_p_output)
  
  beta_p_output <-
    matrix(beta_p_output[, , 1], nrow = nchain, ncol = niter)
  
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

tracePlot_DetectionCovariate <- function(modelResults, index_cov) {
  beta_p_output <- modelResults$modelOutput$beta_p_output
  
  if ((1 + index_cov) <= dim(beta_p_output)[3]) {
    nchain <- nrow(beta_p_output)
    niter <- ncol(beta_p_output)
    
    beta_p_output <-
      matrix(beta_p_output[, , 1 + index_cov], nrow = nchain, ncol = niter)
    
    beta_psi_output_long <- reshape2::melt(beta_p_output)
    
    diagnosticsPlot <-
      ggplot2::ggplot(data = beta_psi_output_long, aes(
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