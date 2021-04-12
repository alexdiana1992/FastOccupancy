#' Simulate data function
#'
#' @param Y number of years 
#' @param S_years number of sites per year
#' @param V_lambda 
#' @param mu_psi 
#' @param b_t 
#' @param beta_psi_tsp 
#' @param beta_psi 
#' @param X 
#' @param a_s_site 
#' @param mu_p 
#' @param beta_p 
#'
#' @return
#' @export
#'
#' @examples
#' simulateData(Y = 3, S_years = rep(50, 3),
#'              V_lambda = 2, mu_psi = 0, b_t = c(1,2,3), beta_psi_tsp = c(0,0),
#'              beta_psi = c(1), X = matrix(runif(50 * 2), 50, 2), a_s_site = rnorm(50),
#'               mu_p = -1, beta_p = 0)
simulateData <- function(Y, S_years, V_lambda, 
                         mu_psi, b_t, beta_psi_tsp, beta_psi, 
                         X, a_s_site,
                         mu_p, beta_p){
  
  # mu_psi <- mu_psi_true
  # b_t <- b_t_true
  # beta_psi_tsp <- beta_psi_tsp_true
  # beta_psi <- beta_psi_true
  # a_s_site <- a_s_site_true
  # mu_p <- mu_p_true
  # beta_p <- beta_p_true
  
  # simulate occupancies
  {
    ncov_psi <- length(beta_psi)
    X_psi_cov <- matrix(stats::rnorm(sum(S_years) * ncov_psi, sd = 1), nrow = sum(S_years), ncol = ncov_psi)
    
    # index_year <- rep(1:Y, each = S)
    index_year <- rep(1:Y, S_years)
    X_psi_year <- data.frame(Year = factor(index_year))
    X_psi_year <- stats::model.matrix(~ . - 1, X_psi_year)
    X_y_index <- apply(X_psi_year, 1, function(x) {which(x != 0)})
    X_year_cov <- scale(as.numeric(X_y_index))
    
    beta_psi_all <- c(mu_psi, b_t, beta_psi_tsp, beta_psi)
    
    # k_s <- rep(1:S, times = Y)
    k_s <- rep(NA, sum(S_years))
    for (y in 1:Y) {
      k_s[sum(S_years[seq_len(y-1)]) + seq_len(S_years[y])] <-  seq_len(S_years[y])
    }
    
    X_psi_tsp_cov <- cbind(X[k_s,1] * X_year_cov, X[k_s,2] * X_year_cov)
    
    X_psi <- cbind(1, X_psi_year, X_psi_tsp_cov, X_psi_cov)
    
    a_s <- a_s_site[k_s]
    
    psi <- logit(X_psi %*% beta_psi_all + a_s)
    
    z <- rep(NA, sum(S_years))
    for (y in 1:Y) {
      for (s in 1:S_years[y]) {
        z[sum(S_years[seq_len(y-1)]) + s] <- stats::rbinom(1, 1, prob = psi[sum(S_years[seq_len(y-1)]) + s])  
      }
    }
    
    visitPerSiteYear <- stats::rpois(sum(S_years), V_lambda)
    z_all <- unlist(lapply(seq_along(z), function(i){
      rep(z[i], visitPerSiteYear[i])
    }))
    # z_all <- lapply(seq_along(z), function(i){
    #   rep(z[i], visitPerSiteYear[i])
    # })
    # z_all <- Reduce("c", z_all)
  }
  
  # simulate detections
  {
    ncov_p <- length(beta_p)
    X_p <- matrix(stats::rnorm(ncov_p * length(z_all)), nrow = length(z_all), ncol = ncov_p)
    
    probsDetection <- logit(mu_p + X_p %*% beta_p)
    
    y_ysv <- rep(NA, length(z_all))
    l <- 1
    for (l in 1:length(z_all)) {
      if(z_all[l] == 0){
        y_ysv[l] <- 0
      } else {
        y_ysv[l] <- stats::rbinom(1, 1, probsDetection[l])
      }
    }
    
  }
  
  # put data in data.frame
  {
    index_site <- rep(NA, sum(S_years))
    for (y in 1:Y) {
      index_site[(sum(S_years[seq_len(y-1)]) + 1:S_years[y])] <- 1:S_years[y]
    }
    
    indexes_site <- unlist(lapply(seq_along(index_site), function(i){
      rep(index_site[i], visitPerSiteYear[i])
    }))
    
    
    indexes_year <- unlist(lapply(seq_along(index_year), function(i){
      rep(index_year[i], visitPerSiteYear[i])
    }))
    
    indexes_covs <- unlist(lapply(seq_len(nrow(X_psi_cov)), function(i){
      rep(i, visitPerSiteYear[i])
    }))
    
    data <- data.frame(
      # Year = rep(index_year, each = V),
      Year = indexes_year,
      Site = indexes_site,
      X = X[indexes_site,],
      # psi_covs = X_psi_cov[rep(1:nrow(X_psi_cov), each = V),],
      psi_covs = X_psi_cov[indexes_covs,],
      p_covs = X_p,
      Occ = y_ysv)
    
  }
  
  return(data)
}

logit <- function(x){
  1 / (1 + exp(-x))
}