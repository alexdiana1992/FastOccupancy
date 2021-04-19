ExtractCovariatesFromText <- function(covariates_text) {
  
  covariates <- sapply(strsplit(covariates_text, split = ",")[[1]], function(x){
    if(grepl("-",x)){
      as.numeric(seq(strsplit(x, split = "-")[[1]][1], strsplit(x, split = "-")[[1]][2]))
    } else {
      as.numeric(x)
    }
  })
  covariates <- as.vector(unlist(covariates))
  if(covariates[1] == 0) covariates <- c()
  
  covariates
}

# true or false if the covariates is categorical or numerical
ExtractClassCovariates <- function(ncov, fac_covariates, column_covariate){
  classCovariates <- rep(T, ncov)
  classCovariates[match(fac_covariates, column_covariate)] <- F
  
  classCovariates
}

Extract_IndexesCovariates <- function(X, column_covariate, ncov, classCovariates) {
  
  indexes_covariates <- c()
  indexes_covariates[1] <- 1
  if(any(column_covariate != 0)){
    
    indexes_covariates <- c()
    indexes_covariates[1] <- 1
    k <- 2
    for (i in 1:ncov) {
      if(classCovariates[i]){
        indexes_covariates[k] <- i + 1
        k <- k + 1
      } else {
        num_levels <- length(unique(X[,i]))
        indexes_covariates[k + 0:(num_levels-2)] <- i + 1
        k <- k + num_levels - 1
      }
    }
    
  }
  
  indexes_covariates
}

logit <- function(x){
  1 / (1 + exp(-x))
}

invLogit <- function(x){
  log(x / (1 - x))
}