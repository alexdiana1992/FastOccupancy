# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

sample_z_cpp <- function(psi, p, k_s) {
    .Call('_FastOccupancy_sample_z_cpp', PACKAGE = 'FastOccupancy', psi, p, k_s)
}

simulateDetections <- function(p, z_all) {
    .Call('_FastOccupancy_simulateDetections', PACKAGE = 'FastOccupancy', p, z_all)
}

computeYearEffect <- function(Y, a_s_unique, beta_psi) {
    .Call('_FastOccupancy_computeYearEffect', PACKAGE = 'FastOccupancy', Y, a_s_unique, beta_psi)
}

sampleBetaPsi <- function(beta, a_s, X, Xbeta, b, invB, n, k, Y, X_centers, ncov_psi, numTimeSpaceCov, X_y_index, X_s_index) {
    .Call('_FastOccupancy_sampleBetaPsi', PACKAGE = 'FastOccupancy', beta, a_s, X, Xbeta, b, invB, n, k, Y, X_centers, ncov_psi, numTimeSpaceCov, X_y_index, X_s_index)
}

sampleBetaP <- function(beta, X, Xbeta, b, invB, ncov_p, p_intercepts, X_y_index, n, k) {
    .Call('_FastOccupancy_sampleBetaP', PACKAGE = 'FastOccupancy', beta, X, Xbeta, b, invB, ncov_p, p_intercepts, X_y_index, n, k)
}

K <- function(x1, x2, a, l) {
    .Call('_FastOccupancy_K', PACKAGE = 'FastOccupancy', x1, x2, a, l)
}

K2 <- function(x1, x2, a, l) {
    .Call('_FastOccupancy_K2', PACKAGE = 'FastOccupancy', x1, x2, a, l)
}

sample_eps_cpp <- function(k_s, sites, c_psi, k, z, Omega, sigma_a) {
    .Call('_FastOccupancy_sample_eps_cpp', PACKAGE = 'FastOccupancy', k_s, sites, c_psi, k, z, Omega, sigma_a)
}

isPointInBandRight <- function(X_tilde, x_grid, y_grid, i, j) {
    .Call('_FastOccupancy_isPointInBandRight', PACKAGE = 'FastOccupancy', X_tilde, x_grid, y_grid, i, j)
}

isPointInBandLeft <- function(X_tilde, x_grid, y_grid, i, j) {
    .Call('_FastOccupancy_isPointInBandLeft', PACKAGE = 'FastOccupancy', X_tilde, x_grid, y_grid, i, j)
}

isPointInBandUp <- function(X_tilde, x_grid, y_grid, i, j) {
    .Call('_FastOccupancy_isPointInBandUp', PACKAGE = 'FastOccupancy', X_tilde, x_grid, y_grid, i, j)
}

isPointInBandDown <- function(X_tilde, x_grid, y_grid, i, j) {
    .Call('_FastOccupancy_isPointInBandDown', PACKAGE = 'FastOccupancy', X_tilde, x_grid, y_grid, i, j)
}

findClosestPoint <- function(XY_sp, X_tilde) {
    .Call('_FastOccupancy_findClosestPoint', PACKAGE = 'FastOccupancy', XY_sp, X_tilde)
}

matrixProductXtOmegaX_year <- function(Y, Omega, X_y_index) {
    .Call('_FastOccupancy_matrixProductXtOmegaX_year', PACKAGE = 'FastOccupancy', Y, Omega, X_y_index)
}

matrixProductXtOmegaX_spatial <- function(X_centers, Omega, X_s_index) {
    .Call('_FastOccupancy_matrixProductXtOmegaX_spatial', PACKAGE = 'FastOccupancy', X_centers, Omega, X_s_index)
}

XpsiYz <- function(X_y_index, z, Y) {
    .Call('_FastOccupancy_XpsiYz', PACKAGE = 'FastOccupancy', X_y_index, z, Y)
}

XpsinoYbetaz <- function(X_s_index, X_centers, ncov, X_cov, beta) {
    .Call('_FastOccupancy_XpsinoYbetaz', PACKAGE = 'FastOccupancy', X_s_index, X_centers, ncov, X_cov, beta)
}

rcpp_log_dmvnorm_fast_cpp <- function(inv_S, diag_S, sigma_s, x) {
    .Call('_FastOccupancy_rcpp_log_dmvnorm_fast_cpp', PACKAGE = 'FastOccupancy', inv_S, diag_S, sigma_s, x)
}

sample_l_grid_cpp <- function(l_s_grid, sigma_s, inv_K_s_grid, diag_K_s_grid, a_l_S, b_l_S, a_s) {
    .Call('_FastOccupancy_sample_l_grid_cpp', PACKAGE = 'FastOccupancy', l_s_grid, sigma_s, inv_K_s_grid, diag_K_s_grid, a_l_S, b_l_S, a_s)
}

