
buildGrid <- function(XY_sp, gridStep){
  
  x_grid <- seq(min(XY_sp[,1]) - (1.5) * gridStep, 
                max(XY_sp[,1]) + (1.5) * gridStep, by = gridStep)
  y_grid <- seq(min(XY_sp[,2]) - (1.5) * gridStep, 
                max(XY_sp[,2]) + (1.5) * gridStep, by = gridStep)
  
  pointInGrid <- matrix(T, nrow = length(x_grid), ncol = length(y_grid))
  
  for (i in 2:(length(x_grid) - 1)) {
    
    for (j in 2:(length(y_grid) - 1)) {
      
      isAnyPointInBandRight <- isPointInBandRight(XY_sp, x_grid, y_grid, i - 1, j - 1)
      
      isAnyPointInBandLeft <- isPointInBandLeft(XY_sp, x_grid, y_grid, i - 1, j - 1)
      
      isAnyPointInBandUp <- isPointInBandUp(XY_sp, x_grid, y_grid, i - 1, j - 1)
      
      isAnyPointInBandDown <- isPointInBandDown(XY_sp, x_grid, y_grid, i - 1, j - 1)
      
      if(!isAnyPointInBandRight | !isAnyPointInBandLeft | !isAnyPointInBandUp | !isAnyPointInBandDown){
        pointInGrid[i,j] <- F
      }
      
    } 
    
  }
  
  pointInGrid <- pointInGrid[-c(1,nrow(pointInGrid)),]
  pointInGrid <- pointInGrid[,-c(1,ncol(pointInGrid))]
  x_grid <- x_grid[-c(1,length(x_grid))]
  y_grid <- y_grid[-c(1,length(y_grid))]
  
  allPoints <- cbind(expand.grid(x_grid, y_grid), as.vector((pointInGrid)))
  allPoints <- allPoints[allPoints[,3],-3]  
  
  allPoints
}


#' Build spatial grid
#' 
#' @description This function shows the spatial grid used in the approximation of the site-specific autocorrelated random effects for a particular value of the grid step.
#' @param X_cov The variable containing the x coordinates of the sites.
#' @param Y_cov The variable containing the y coordinates of the sites.
#' @param gridStep Value of the grid step.
#' 
#' @importFrom magrittr %>%
#' 
#' @export
#' 
#' @return Plot of the grid with number of points used in the approximation.
#' 
#' @examples
#' 
#' buildSpatialGrid(sampleData$X, sampleData$Y, gridStep = .3)
#' 
buildSpatialGrid <- function(X_cov, Y_cov, gridStep){
  
  meanXsp <- mean(X_cov)
  sdXsp <- stats::sd(X_cov)
  meanYsp <- mean(Y_cov)
  sdYsp <- stats::sd(Y_cov)
  X_sp <- (X_cov - mean(X_cov)) / stats::sd(X_cov)
  Y_sp <- (Y_cov - mean(Y_cov)) / stats::sd(Y_cov)
  
  # cluster the sites
  
  uniqueIndexesSite <- which(!duplicated(cbind(X_cov,Y_cov)))
  X_sp_unique <- X_sp[uniqueIndexesSite]
  Y_sp_unique <- Y_sp[uniqueIndexesSite]
  
  # numCenters <- 400
  # kmeansclustering <- kmeans(cbind(X_sp_unique, Y_sp_unique), centers = numCenters, iter.max = 100)
  # X_tilde <- kmeansclustering$centers
  # X_psi$Center <- as.factor(kmeansclustering$cluster[k_s$Site])
  # X_psi_s <- model.matrix( ~ . - 1, data = X_psi[,c("Center"),drop = F])
  # X_s_index <- kmeansclustering$cluster[k_s$Site]
  
  XY_sp <- cbind(X_sp_unique, Y_sp_unique)
  X_tilde <- as.matrix(buildGrid(XY_sp, gridStep))
  XY_centers <- findClosestPoint(XY_sp, X_tilde)
  X_tilde <- X_tilde[unique(XY_centers),]
  XY_centers <- findClosestPoint(XY_sp, X_tilde)
  numCenters <- nrow(X_tilde)
  
  X_tilde_original <- X_tilde
  X_tilde_original[,1] <- X_tilde[,1] * sdXsp + meanXsp
  X_tilde_original[,2] <- X_tilde[,2] * sdXsp + meanXsp
  
  ggplot2::ggplot() + 
    ggplot2::geom_point(data = NULL, ggplot2::aes(x = X_tilde_original[,1],
                                                  y = X_tilde_original[,2]), color = "black", size = 3, shape = 15) +
    ggplot2::ylab("Y") + ggplot2::scale_x_continuous(name = "X") + 
    ggplot2::ggtitle(label = paste("Number of points = ",nrow(X_tilde))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 13),
                   axis.title = ggplot2::element_text(size = 13, face = "bold"),
                   axis.text = ggplot2::element_text(size = 13, face = "bold", angle = 90),
                   panel.grid.major = ggplot2::element_line(colour="grey", size=0.015),
                   panel.background = ggplot2::element_rect(fill = "white", color = "black")) 
  
}