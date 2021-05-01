
update_p <- function(beta_p, y, z_all, X_p, b_p, inv_B_p, Xbeta_p, 
                     ncov_p, p_intercepts, X_y_index_p, usingYearDetProb){
  
  k <- y - .5
  n <- rep(1, length(k))
  
  k <- k[z_all==1]
  n <- n[z_all==1]
  X_p_present <- X_p[z_all == 1,,drop = FALSE]
  X_y_index_p_present <- X_y_index_p[z_all == 1]

  Xbeta <- Xbeta_p[z_all == 1]
  
  # list_beta_p <- sample_beta_omega_cpp_parallel(beta_p, X_p_present, b_p, B_p, n, k)
  # list_beta_p <- sample_beta_omega_cpp(beta_p, X_p_present, b_p, inv_B_p, n, k)
  # beta_p <- list_beta_p$beta
  beta_p <- sampleBetaP(beta_p, X_p_present, Xbeta, b_p, inv_B_p,
                          ncov_p, p_intercepts, X_y_index_p_present, n, k)
  
  Xbeta <- X_p[,-(1:p_intercepts),drop = F] %*% beta_p[-(1:p_intercepts)]
  XbetaY <- beta_p[X_y_index_p]
  Xbeta <- Xbeta + XbetaY

  p <- as.vector(logit(Xbeta))
  
  # p <- as.vector(logit(X_p %*% beta_p))
  
  list("p" = p, "beta" = beta_p, "Xbeta" = Xbeta)
}