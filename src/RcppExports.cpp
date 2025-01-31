// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sample_z_cpp
arma::vec sample_z_cpp(arma::vec psi, arma::vec p, arma::mat k_s);
RcppExport SEXP _FastOccupancy_sample_z_cpp(SEXP psiSEXP, SEXP pSEXP, SEXP k_sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type k_s(k_sSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_z_cpp(psi, p, k_s));
    return rcpp_result_gen;
END_RCPP
}
// simulateDetections
arma::vec simulateDetections(arma::vec p, arma::vec z_all);
RcppExport SEXP _FastOccupancy_simulateDetections(SEXP pSEXP, SEXP z_allSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z_all(z_allSEXP);
    rcpp_result_gen = Rcpp::wrap(simulateDetections(p, z_all));
    return rcpp_result_gen;
END_RCPP
}
// computeYearEffect
arma::vec computeYearEffect(int Y, arma::vec a_s_unique, arma::vec beta_psi);
RcppExport SEXP _FastOccupancy_computeYearEffect(SEXP YSEXP, SEXP a_s_uniqueSEXP, SEXP beta_psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a_s_unique(a_s_uniqueSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type beta_psi(beta_psiSEXP);
    rcpp_result_gen = Rcpp::wrap(computeYearEffect(Y, a_s_unique, beta_psi));
    return rcpp_result_gen;
END_RCPP
}
// sampleBetaPsi
List sampleBetaPsi(arma::vec beta, arma::vec a_s, arma::mat& X, arma::vec Xbeta, arma::vec b, arma::mat invB, arma::vec n, arma::vec k, int Y, int X_centers, int ncov_psi, int numTimeSpaceCov, IntegerVector X_y_index, IntegerVector X_s_index);
RcppExport SEXP _FastOccupancy_sampleBetaPsi(SEXP betaSEXP, SEXP a_sSEXP, SEXP XSEXP, SEXP XbetaSEXP, SEXP bSEXP, SEXP invBSEXP, SEXP nSEXP, SEXP kSEXP, SEXP YSEXP, SEXP X_centersSEXP, SEXP ncov_psiSEXP, SEXP numTimeSpaceCovSEXP, SEXP X_y_indexSEXP, SEXP X_s_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a_s(a_sSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xbeta(XbetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invB(invBSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type X_centers(X_centersSEXP);
    Rcpp::traits::input_parameter< int >::type ncov_psi(ncov_psiSEXP);
    Rcpp::traits::input_parameter< int >::type numTimeSpaceCov(numTimeSpaceCovSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type X_y_index(X_y_indexSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type X_s_index(X_s_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleBetaPsi(beta, a_s, X, Xbeta, b, invB, n, k, Y, X_centers, ncov_psi, numTimeSpaceCov, X_y_index, X_s_index));
    return rcpp_result_gen;
END_RCPP
}
// sampleBetaP
arma::vec sampleBetaP(arma::vec beta, arma::mat& X, arma::vec Xbeta, arma::vec b, arma::mat invB, int ncov_p, int p_intercepts, IntegerVector X_y_index, arma::vec n, arma::vec k);
RcppExport SEXP _FastOccupancy_sampleBetaP(SEXP betaSEXP, SEXP XSEXP, SEXP XbetaSEXP, SEXP bSEXP, SEXP invBSEXP, SEXP ncov_pSEXP, SEXP p_interceptsSEXP, SEXP X_y_indexSEXP, SEXP nSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Xbeta(XbetaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type invB(invBSEXP);
    Rcpp::traits::input_parameter< int >::type ncov_p(ncov_pSEXP);
    Rcpp::traits::input_parameter< int >::type p_intercepts(p_interceptsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type X_y_index(X_y_indexSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleBetaP(beta, X, Xbeta, b, invB, ncov_p, p_intercepts, X_y_index, n, k));
    return rcpp_result_gen;
END_RCPP
}
// K
arma::mat K(arma::vec x1, arma::vec x2, double a, double l);
RcppExport SEXP _FastOccupancy_K(SEXP x1SEXP, SEXP x2SEXP, SEXP aSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(K(x1, x2, a, l));
    return rcpp_result_gen;
END_RCPP
}
// K2
arma::mat K2(arma::mat x1, arma::mat x2, double a, double l);
RcppExport SEXP _FastOccupancy_K2(SEXP x1SEXP, SEXP x2SEXP, SEXP aSEXP, SEXP lSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x2(x2SEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type l(lSEXP);
    rcpp_result_gen = Rcpp::wrap(K2(x1, x2, a, l));
    return rcpp_result_gen;
END_RCPP
}
// sample_eps_cpp
arma::vec sample_eps_cpp(arma::vec k_s, arma::vec sites, arma::vec& c_psi, arma::vec k, arma::vec z, arma::vec Omega, double sigma_a);
RcppExport SEXP _FastOccupancy_sample_eps_cpp(SEXP k_sSEXP, SEXP sitesSEXP, SEXP c_psiSEXP, SEXP kSEXP, SEXP zSEXP, SEXP OmegaSEXP, SEXP sigma_aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type k_s(k_sSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sites(sitesSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type c_psi(c_psiSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_a(sigma_aSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_eps_cpp(k_s, sites, c_psi, k, z, Omega, sigma_a));
    return rcpp_result_gen;
END_RCPP
}
// isPointInBandRight
bool isPointInBandRight(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j);
RcppExport SEXP _FastOccupancy_isPointInBandRight(SEXP X_tildeSEXP, SEXP x_gridSEXP, SEXP y_gridSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_tilde(X_tildeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_grid(y_gridSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(isPointInBandRight(X_tilde, x_grid, y_grid, i, j));
    return rcpp_result_gen;
END_RCPP
}
// isPointInBandLeft
bool isPointInBandLeft(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j);
RcppExport SEXP _FastOccupancy_isPointInBandLeft(SEXP X_tildeSEXP, SEXP x_gridSEXP, SEXP y_gridSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_tilde(X_tildeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_grid(y_gridSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(isPointInBandLeft(X_tilde, x_grid, y_grid, i, j));
    return rcpp_result_gen;
END_RCPP
}
// isPointInBandUp
bool isPointInBandUp(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j);
RcppExport SEXP _FastOccupancy_isPointInBandUp(SEXP X_tildeSEXP, SEXP x_gridSEXP, SEXP y_gridSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_tilde(X_tildeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_grid(y_gridSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(isPointInBandUp(X_tilde, x_grid, y_grid, i, j));
    return rcpp_result_gen;
END_RCPP
}
// isPointInBandDown
bool isPointInBandDown(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j);
RcppExport SEXP _FastOccupancy_isPointInBandDown(SEXP X_tildeSEXP, SEXP x_gridSEXP, SEXP y_gridSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X_tilde(X_tildeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type x_grid(x_gridSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y_grid(y_gridSEXP);
    Rcpp::traits::input_parameter< int >::type i(iSEXP);
    Rcpp::traits::input_parameter< int >::type j(jSEXP);
    rcpp_result_gen = Rcpp::wrap(isPointInBandDown(X_tilde, x_grid, y_grid, i, j));
    return rcpp_result_gen;
END_RCPP
}
// findClosestPoint
IntegerVector findClosestPoint(arma::mat XY_sp, arma::mat X_tilde);
RcppExport SEXP _FastOccupancy_findClosestPoint(SEXP XY_spSEXP, SEXP X_tildeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type XY_sp(XY_spSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type X_tilde(X_tildeSEXP);
    rcpp_result_gen = Rcpp::wrap(findClosestPoint(XY_sp, X_tilde));
    return rcpp_result_gen;
END_RCPP
}
// matrixProductXtOmegaX_year
arma::mat matrixProductXtOmegaX_year(int Y, arma::vec Omega, arma::vec X_y_index);
RcppExport SEXP _FastOccupancy_matrixProductXtOmegaX_year(SEXP YSEXP, SEXP OmegaSEXP, SEXP X_y_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type X_y_index(X_y_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixProductXtOmegaX_year(Y, Omega, X_y_index));
    return rcpp_result_gen;
END_RCPP
}
// matrixProductXtOmegaX_spatial
arma::mat matrixProductXtOmegaX_spatial(int X_centers, arma::vec Omega, arma::vec X_s_index);
RcppExport SEXP _FastOccupancy_matrixProductXtOmegaX_spatial(SEXP X_centersSEXP, SEXP OmegaSEXP, SEXP X_s_indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type X_centers(X_centersSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type Omega(OmegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type X_s_index(X_s_indexSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixProductXtOmegaX_spatial(X_centers, Omega, X_s_index));
    return rcpp_result_gen;
END_RCPP
}
// XpsiYz
arma::vec XpsiYz(arma::vec X_y_index, arma::vec z, int Y);
RcppExport SEXP _FastOccupancy_XpsiYz(SEXP X_y_indexSEXP, SEXP zSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X_y_index(X_y_indexSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(XpsiYz(X_y_index, z, Y));
    return rcpp_result_gen;
END_RCPP
}
// XpsinoYbetaz
arma::vec XpsinoYbetaz(arma::vec X_s_index, int X_centers, int ncov, arma::mat& X_cov, arma::vec& beta);
RcppExport SEXP _FastOccupancy_XpsinoYbetaz(SEXP X_s_indexSEXP, SEXP X_centersSEXP, SEXP ncovSEXP, SEXP X_covSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type X_s_index(X_s_indexSEXP);
    Rcpp::traits::input_parameter< int >::type X_centers(X_centersSEXP);
    Rcpp::traits::input_parameter< int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type X_cov(X_covSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(XpsinoYbetaz(X_s_index, X_centers, ncov, X_cov, beta));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_log_dmvnorm_fast_cpp
double rcpp_log_dmvnorm_fast_cpp(arma::mat& inv_S, arma::vec& diag_S, double sigma_s, arma::vec& x);
RcppExport SEXP _FastOccupancy_rcpp_log_dmvnorm_fast_cpp(SEXP inv_SSEXP, SEXP diag_SSEXP, SEXP sigma_sSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type inv_S(inv_SSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type diag_S(diag_SSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_s(sigma_sSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpp_log_dmvnorm_fast_cpp(inv_S, diag_S, sigma_s, x));
    return rcpp_result_gen;
END_RCPP
}
// sample_l_grid_cpp
arma::vec sample_l_grid_cpp(arma::vec l_s_grid, double sigma_s, arma::cube& inv_K_s_grid, arma::mat& diag_K_s_grid, double a_l_S, double b_l_S, arma::vec a_s);
RcppExport SEXP _FastOccupancy_sample_l_grid_cpp(SEXP l_s_gridSEXP, SEXP sigma_sSEXP, SEXP inv_K_s_gridSEXP, SEXP diag_K_s_gridSEXP, SEXP a_l_SSEXP, SEXP b_l_SSEXP, SEXP a_sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type l_s_grid(l_s_gridSEXP);
    Rcpp::traits::input_parameter< double >::type sigma_s(sigma_sSEXP);
    Rcpp::traits::input_parameter< arma::cube& >::type inv_K_s_grid(inv_K_s_gridSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type diag_K_s_grid(diag_K_s_gridSEXP);
    Rcpp::traits::input_parameter< double >::type a_l_S(a_l_SSEXP);
    Rcpp::traits::input_parameter< double >::type b_l_S(b_l_SSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type a_s(a_sSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_l_grid_cpp(l_s_grid, sigma_s, inv_K_s_grid, diag_K_s_grid, a_l_S, b_l_S, a_s));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FastOccupancy_sample_z_cpp", (DL_FUNC) &_FastOccupancy_sample_z_cpp, 3},
    {"_FastOccupancy_simulateDetections", (DL_FUNC) &_FastOccupancy_simulateDetections, 2},
    {"_FastOccupancy_computeYearEffect", (DL_FUNC) &_FastOccupancy_computeYearEffect, 3},
    {"_FastOccupancy_sampleBetaPsi", (DL_FUNC) &_FastOccupancy_sampleBetaPsi, 14},
    {"_FastOccupancy_sampleBetaP", (DL_FUNC) &_FastOccupancy_sampleBetaP, 10},
    {"_FastOccupancy_K", (DL_FUNC) &_FastOccupancy_K, 4},
    {"_FastOccupancy_K2", (DL_FUNC) &_FastOccupancy_K2, 4},
    {"_FastOccupancy_sample_eps_cpp", (DL_FUNC) &_FastOccupancy_sample_eps_cpp, 7},
    {"_FastOccupancy_isPointInBandRight", (DL_FUNC) &_FastOccupancy_isPointInBandRight, 5},
    {"_FastOccupancy_isPointInBandLeft", (DL_FUNC) &_FastOccupancy_isPointInBandLeft, 5},
    {"_FastOccupancy_isPointInBandUp", (DL_FUNC) &_FastOccupancy_isPointInBandUp, 5},
    {"_FastOccupancy_isPointInBandDown", (DL_FUNC) &_FastOccupancy_isPointInBandDown, 5},
    {"_FastOccupancy_findClosestPoint", (DL_FUNC) &_FastOccupancy_findClosestPoint, 2},
    {"_FastOccupancy_matrixProductXtOmegaX_year", (DL_FUNC) &_FastOccupancy_matrixProductXtOmegaX_year, 3},
    {"_FastOccupancy_matrixProductXtOmegaX_spatial", (DL_FUNC) &_FastOccupancy_matrixProductXtOmegaX_spatial, 3},
    {"_FastOccupancy_XpsiYz", (DL_FUNC) &_FastOccupancy_XpsiYz, 3},
    {"_FastOccupancy_XpsinoYbetaz", (DL_FUNC) &_FastOccupancy_XpsinoYbetaz, 5},
    {"_FastOccupancy_rcpp_log_dmvnorm_fast_cpp", (DL_FUNC) &_FastOccupancy_rcpp_log_dmvnorm_fast_cpp, 4},
    {"_FastOccupancy_sample_l_grid_cpp", (DL_FUNC) &_FastOccupancy_sample_l_grid_cpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_FastOccupancy(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
