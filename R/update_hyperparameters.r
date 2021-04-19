
rinvgamma <- function(a, b){
  1 / stats::rgamma(1, shape = a, rate = b)
}


dinvgamma <- function (x, shape, scale = 1)  {
  log.density <- shape * log(scale) - lgamma(shape) - 
    (shape + 1) * log(x) - (scale /  x)
  return(exp(log.density))
}


update_sigma <- function(a_s, k_s, a_sigma, b_sigma){
  
  a_s_sites <- a_s[!duplicated(k_s$Site)]
  n <- length(a_s_sites)
  
  rinvgamma(a_sigma + (n / 2), b_sigma + sum(a_s_sites^2)/2)
  
}

update_l_sigma_integrated <- function(l_T, sigma_T, a_l_T, b_l_T, a_ig, b_ig, Y, sigma_psi,
                                      z, XbetaY, Xbetas, Xbeta_cov, eps_s, Omega,  
                                      b_psi_b, sd_l, sd_sigma_T,  X_y_index, usingSpatial){
  
  l_T_star <- stats::rnorm(1, l_T, sd_l)
  sigma_T_star <- stats::rnorm(1, sigma_T, sd_sigma_T)
  
  if(l_T_star < 4 & l_T_star > 0 & sigma_T_star > 0){
    
    K_l <- K(1:Y, 1:Y, sigma_T^2, l_T) + sigma_psi^2 + diag(exp(-10), nrow = Y)
    K_l_star <- K(1:Y, 1:Y, sigma_T_star^2, l_T_star) + sigma_psi^2 + diag(exp(-10), nrow = Y)
    
    # term 1
    
    tXOmegX <- matrixProductXtOmegaX_year(Y, Omega, X_y_index)
    
    logterm1_1 <- -.5 * log(FastGP::rcppeigen_get_det(K_l_star %*% tXOmegX + diag(1, nrow = Y)))
    logterm1_2 <- -.5 * log(FastGP::rcppeigen_get_det(K_l %*% tXOmegX + diag(1, nrow = Y)))
    
    logdeterminantsRatio <- logterm1_1 - logterm1_2
    
    # term 2
    
    logmuKlmu <- - .5 * b_psi_b %*% (FastGP::rcppeigen_invert_matrix(K_l_star) - 
                                       FastGP::rcppeigen_invert_matrix(K_l)) %*% b_psi_b
    
    # term 3
    
    if(usingSpatial){
      c_i <- Xbetas + Xbeta_cov + eps_s  
    } else {
      c_i <- Xbeta_cov + eps_s
    }
    
    
    k_i <- z - .5
    z_tilde <- k_i - c_i * as.vector(Omega)
    
    Xtz_tilde <- XpsiYz(X_y_index, z_tilde, Y)
    
    mu_tilde <- Xtz_tilde + FastGP::rcppeigen_invert_matrix(K_l) %*% b_psi_b
    mu_tilde_star <- Xtz_tilde + FastGP::rcppeigen_invert_matrix(K_l_star) %*% b_psi_b
    
    XtOXpB <- tXOmegX + FastGP::rcppeigen_invert_matrix(K_l)
    XtOXpB_star <- tXOmegX + FastGP::rcppeigen_invert_matrix(K_l_star)
    
    log_term2ratio_star <- .5 * ( t(mu_tilde_star) %*% FastGP::rcppeigen_invert_matrix(XtOXpB_star) %*% mu_tilde_star )
    log_term2ratio <- .5 * ( t(mu_tilde) %*% FastGP::rcppeigen_invert_matrix(XtOXpB) %*% mu_tilde )
    
    logterm3 <- log_term2ratio_star - log_term2ratio
    
    likelihood <- exp(logdeterminantsRatio + logmuKlmu +  logterm3)
    
    prior_l <- exp(stats::dgamma(l_T_star, a_l_T, b_l_T, log = T) - stats::dgamma(l_T, a_l_T, b_l_T, log = T))
    prior_sigma <- exp(dinvgamma(sigma_T_star, a_ig, b_ig) - dinvgamma(sigma_T, a_ig, b_ig))
    
    posterior <-  likelihood * prior_l * prior_sigma
    
    if(!is.na(posterior)){
      if(stats::runif(1) < posterior){
        l_T <- l_T_star
        sigma_T <- sigma_T_star
      }  
    }
    
  }
  
  list("l_T" = l_T, "sigma_T" = sigma_T)
}

rcpp_log_dmvnorm_fast <- function (inv_S, diag_S, sigma_s, x) {
  n <- length(x)
  # return((-n/2) * log(2 * pi) - sum(log(diag_S)) - (1/2) * x %*% inv_S %*% x)
  return((-n/2) * log(2 * pi) - sum(log(diag_S) + log(sigma_s)) -
           (1/2) * (1 / sigma_s^2) * x %*% inv_S %*% x)
}

sample_l_grid <- function(l_s_grid, inv_K_s_grid, diag_K_s_grid,
                          a_l_s, b_l_s, a_s, sigma_s){
  
  posterior_val <- rep(NA, length(l_s_grid))
  
  for (j in 1:length(l_s_grid)) {
    
    l_s <- l_s_grid[j]
    # K_s_grid_j <- K_s_grid[,,j] * sigma_s^2
    # inv_K_s_grid_j <- inv_K_s_grid[,,j] / sigma_s^2
    # diag_K_s_grid_j <- diag_K_s_grid[,j] * sigma_s
    
    # loglikelihood <- rcpp_log_dmvnorm_fast(inv_K_s_grid[,,j],
    #                                        diag_K_s_grid[,j], sigma_s, a_s)
    
    loglikelihood <- ((-length(a_s)/2) * log(2 * pi) - sum(log(diag_K_s_grid[,j]) + log(sigma_s)) -
                        (1/2) * (1 / sigma_s^2) * a_s %*% inv_K_s_grid[,,j] %*% a_s)
    
    # loglikelihood <- rcpp_log_dmvnorm_fast(1, inv_K_s_grid_j,
    # diag_K_s_grid_j,  a_s)
    
    # Sigma_l <- K2(X_tilde, X_tilde, sigma_s^2, l_s) + diag(exp(-10), nrow = nrow(X_tilde))
    
    # (loglikelihood2 <- rcpp_log_dmvnorm( Sigma_l, rep(0, X_centers), a_s, F))
    
    logPrior <- dgamma(l_s, a_l_s, b_l_s, log = T)
    
    posterior_val[j] <- logPrior + loglikelihood
    
  }
  
  # posterior_val <- sample_l_grid_cpp(l_s_grid, sigma_s, inv_K_s_grid, diag_K_s_grid,
  # a_l_s, b_l_s, a_s)
  
  posterior_val <- posterior_val - max(posterior_val)
  
  index_l_grid <- sample(1:length(l_s_grid), size = 1, 
                         prob = exp(posterior_val) / sum(exp(posterior_val)))
  
  index_l_grid
}

update_hyperparameters <- function(l_T, a_l_T, b_l_T, sd_l_T, sd_sigma_T,
                                   sigma_T, a_sigma_T, b_sigma_T, Y,
                                   beta_psi, inv_B_psi, 
                                   b_psi, sigma_psi,
                                   l_s_grid, K_s_grid, inv_K_s_grid, 
                                   inv_chol_K_s_grid, diag_K_s_grid,
                                   a_l_s, b_l_s, 
                                   sigma_s, a_sigma_s, b_sigma_s, X_tilde,
                                   a_sigma_eps, b_sigma_eps,
                                   usingSpatial, 
                                   XbetaY, Xbetas, Xbeta_cov,
                                   eps_s, k_s, 
                                   z, X_psi, Omega, X_y_index){
  
  
  # update l_T and sigma_T -------------------
  
  list_l_sigma <- update_l_sigma_integrated(l_T, sigma_T, a_l_T, b_l_T, a_sigma_T, b_sigma_T, Y, 
                                            sigma_psi,
                                            z, betaY, Xbetas, Xbeta_cov, eps_s, Omega,  
                                            b_psi[1:Y], sd_l = sd_l_T, 
                                            sd_sigma_T, X_y_index, usingSpatial)
  l_T <- list_l_sigma$l_T
  sigma_T <- list_l_sigma$sigma_T
  
  # reupdate inv_B_psi ------------
  
  K_l <- K(1:Y, 1:Y, sigma_T^2, l_T) + sigma_psi^2 + diag(exp(-8), nrow = Y)
  
  inverse_Kl <- FastGP::rcppeigen_invert_matrix(K_l)
  inverse_Kl[lower.tri(K_l)] <- t(inverse_Kl)[lower.tri(K_l)]
  
  inv_B_psi[1:Y, 1:Y] <- inverse_Kl
  
  # update l_s -------------------
  
  if(usingSpatial){
    index_ls_grid <- sample_l_grid(l_s_grid, inv_K_s_grid, diag_K_s_grid,
                                   a_l_s, b_l_s, beta_psi[Y + 1:nrow(X_tilde)], sigma_s)
    l_s <- l_s_grid[index_ls_grid]
  } else {
    l_s <- 0
    index_ls_grid <- 0
  }
  
  # update sigma_s ---------------------
  
  if(usingSpatial){
    
    X_centers <- nrow(X_tilde)
    
    inv_chol_Kl <- inv_chol_K_s_grid[,,index_ls_grid]
    
    a_s <- beta_psi[Y + 1:X_centers]
    
    a_ig <- X_centers / 2
    Ltas <- inv_chol_Kl %*% a_s
    b_ig <- (t(Ltas) %*% Ltas) / 2
    
    sigma_s <- sqrt(rinvgamma(a_sigma_s + a_ig, b_sigma_s + b_ig))
  } 
  
  # reupdate inv_B_psi ------------
  
  if(usingSpatial){
    X_centers <- nrow(X_tilde)
    
    inv_K_lsigma <- inv_K_s_grid[,,index_ls_grid] / sigma_s^2
    
    inv_B_psi[Y + 1:X_centers, Y + 1:X_centers] <- inv_K_lsigma  
  }
  
  # update sigma_as ------------------------------
  
  sigma_eps <- sqrt(update_sigma(eps_s, k_s, a_sigma_eps, b_sigma_eps))
  
  list("l_T" = l_T,
       "sigma_T" = sigma_T,
       "inv_B_psi" = inv_B_psi,
       "l_s" = l_s,
       "index_ls_grid" = index_ls_grid,
       "sigma_s" = sigma_s,
       "sigma_eps" = sigma_eps)
}
