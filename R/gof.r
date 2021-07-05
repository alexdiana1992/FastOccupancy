

#' Goodness of fit of the yearly detection
#' 
#' @description This function performs a goodness of fit of the yearly detection
#' as described in Diana et al. (2021).
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A plot with the goddness of fit. The red points represent the true counts and
#' the interval represents the posterior of the yearly counts as estimated from the model 
#' 
#' @examples
#' 
#' plotGOFYearDetections(sampleResults)
#' 
plotGOFYearDetections <- function(modelResults) {
  years <- modelResults$dataCharacteristics$Years
  Y <- length(years)
  
  gofYear_output <-
    modelResults$modelOutput$GOF_output$gofYear_output
  gofYear_output <- apply(gofYear_output, 3, c)
  
  CI_gofs <- apply(gofYear_output, 2, function(x) {
    c(stats::quantile(x, probs = c(0.025, 0.975)), mean(x))
    
  })
  
  trueYearStatistics <-
    modelResults$modelOutput$GOF_output$trueYearStatistics
  
  
  ggplot2::ggplot(
    data = NULL,
    ggplot2::aes(
      x = years,
      y = trueYearStatistics$Detections,
      ymin = CI_gofs[1,],
      ymax = CI_gofs[2,]
    )
  ) +
    ggplot2::geom_errorbar() + ggplot2::geom_line() + ggplot2::geom_point(color = "red") +
    ggplot2::ylab("Detections") + ggplot2::scale_x_continuous(name = "Year", breaks = years) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
      axis.title = ggplot2::element_text(size = 20, face = "bold"),
      axis.text = ggplot2::element_text(
        size = 13,
        face = "bold",
        angle = 90
      ),
      panel.grid.major = ggplot2::element_line(colour = "grey", size =
                                                 0.15),
      panel.background = ggplot2::element_rect(fill = "white", color = "black")
    )
  
  
}

#' Goodness of fit of the detection in each region
#' 
#' @description This function performs a goodness of fit of the detection
#' in each region as described in Diana et al. (2021).
#' @param modelResults Output of the function \code{runModel}
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return A map with the spatial points used in the approximation. The color
#' used varies according to wether the true value of statistcs is inside, above 
#' or below the 95\% credible interval.
#' 
#' @examples
#' 
#' plotGOFSpaceDections(sampleResults)
#' 
plotGOFSpaceDections <- function(modelResults) {
  
  usingSpatial <- modelResults$dataCharacteristics$usingSpatial
  
  if(usingSpatial){
    
    X_tilde <- modelResults$dataCharacteristics$X_tilde
    X_centers <- nrow(X_tilde)
    
    years <- modelResults$dataCharacteristics$Years
    Y <- length(years)
    
    gofSpace_output <-
      modelResults$modelOutput$GOF_output$gofSpace_output
    gofSpace_output <- apply(gofSpace_output, 3, c)
    
    CI_gofs <- apply(gofSpace_output, 2, function(x) {
      stats::quantile(x, probs = c(0.025, 0.975))
    })
    
    trueSpaceStatistics <-
      modelResults$modelOutput$GOF_output$trueSpaceStatistics
    
    GOF_interval <-
      ifelse(
        trueSpaceStatistics$Detections < CI_gofs[1, ],
        "OVER",
        ifelse(trueSpaceStatistics$Detections > CI_gofs[2, ], "UNDER", "INSIDE")
      )
    
    ggplot2::ggplot(data = NULL, ggplot2::aes(x = X_tilde[, 1],
                                              y = X_tilde[, 2],
                                              color = GOF_interval)) +  ggplot2::geom_point(shape = 15, size = 2) +
      ggplot2::ylab("Y") + ggplot2::xlab("X") + ggplot2::guides(color = ggplot2::guide_legend(title = "GOF")) +
      ggplot2::scale_color_manual(values = c("#87D146", "red", "orange")) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 20),
        axis.title = ggplot2::element_text(size = 20, face = "bold"),
        axis.text = ggplot2::element_text(
          size = 13,
          face = "bold",
          angle = 90
        ),
        panel.grid.major = ggplot2::element_line(colour = "grey", size = 0.15),
        panel.background = ggplot2::element_rect(fill = "white", color = "black")
      )
    
  } else {
    
    print("No spatial effect present")
    
  }
  
  
  
  # idx_Sites_order <- order(-trueSpaceStatistics$Detections)
  # idx_Sites <- idx_Sites_order[1:100]
  # # idx_Sites <- 1:100
  # 
  # ggplot(data = NULL, ggplot2::aes(x = (1:nrow(X_tilde))[idx_Sites],
  #                         y = trueSpaceStatistics$Detections[idx_Sites],
  #                         ymin = CI_gofs[1,][idx_Sites],
  #                         ymax = CI_gofs[2,][idx_Sites])) +
  #   geom_errorbar() + geom_point(color = "red", size = 1) +
  #   ylim(c(0, max(trueSpaceStatistics$Detections))) +
  #   ylab("Detections") + scale_x_continuous(name = "Site", breaks = years) +
  #   theme(plot.title = element_text(hjust = 0.5, size = 20),
  #         axis.title = element_text(size = 20, face = "bold"),
  #         axis.text = element_text(size = 13, face = "bold", angle = 90),
  #         panel.grid.major = element_line(colour="grey", size=0.2),
  #         panel.background = element_rect(fill = "white", color = "black"))
  
  
  
}
