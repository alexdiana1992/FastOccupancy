#include "RcppArmadillo.h"

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

const double TRUNC = .64;
const double TRUNC_RECIP = 1.0 / .64;
const double log2pi = std::log(2.0 * M_PI);

// Mathematical constants computed using Wolfram Alpha
#define MATH_PI        3.141592653589793238462643383279502884197169399375105820974
#define MATH_PI_2      1.570796326794896619231321691639751442098584699687552910487
#define MATH_2_PI      0.636619772367581343075535053490057448137838582961825794990
#define MATH_PI2       9.869604401089358618834490999876151135313699407240790626413
#define MATH_PI2_2     4.934802200544679309417245499938075567656849703620395313206
#define MATH_SQRT1_2   0.707106781186547524400844362104849039284835937688474036588
#define MATH_SQRT_PI_2 1.253314137315500251207882642405522626503493370304969158314
#define MATH_LOG_PI    1.144729885849400174143427351353058711647294812915311571513
#define MATH_LOG_2_PI  -0.45158270528945486472619522989488214357179467855505631739
#define MATH_LOG_PI_2  0.451582705289454864726195229894882143571794678555056317392

///////// SAMPLE LATENT OCCUPANCIES

// [[Rcpp::export]]
arma::vec sample_z_cpp(arma::vec psi, arma::vec p, arma::mat k_s){
  
  arma::vec z = arma::zeros(k_s.n_rows);
  
  // this loops over p
  int l = 0;
  
  for(int i = 0; (unsigned)i < z.size(); i++){
    
    if(k_s(i, 0) == 1){
      
      z[i] = 1;
      l += k_s(i,1);
      
    } else {
      
      double prod1mp = 1; //prod(1 - p[i + seq_len(k_s[i,4])])
      for(int k = 0; k < k_s(i,1); k++){
        prod1mp *= (1 - p(l + k));
      }
      
      double p_zsequal1 = (psi[i] * prod1mp) / (psi[i] * prod1mp + (1 - psi[i])) ;
      z[i] = R::rbinom(1, p_zsequal1);
      l += k_s(i,1);
      
    }
    
  }
  
  return(z);
}  

///////// GOODNESS OF FIT

// [[Rcpp::export]]
arma::vec simulateDetections(arma::vec p, arma::vec z_all){
  
  arma::vec y = arma::zeros(p.size());
  
  // this loops over p
  for(int i = 0; (unsigned)i < y.size(); i++){
    
    if(z_all(i) == 1){
      
      y[i] = R::rbinom(1, p[i]);
      
    } 
    
  }
  
  return(y);
}  

arma::vec logit(arma::vec x){
  return(1 / (1 + exp(-x)));  
}

// [[Rcpp::export]]
arma::vec computeYearEffect(int Y, arma::vec a_s_unique, arma::vec beta_psi){
  
  arma::vec yearEffect = arma::zeros(Y);
  
  for(int y = 0; (unsigned)y < Y; y++){
    
    arma::vec occProbs = logit(beta_psi[y] + a_s_unique);
    
    yearEffect[y] = mean(occProbs);
  }
  
  return(yearEffect);
}  

//////// SAMPLE POLYA GAMMA VARIABLES

double aterm(int n, double x, double t) {
  double f = 0;
  if(x <= t) {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
  }
  else {
    f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
  }    
  return (double)exp(f);
}

double exprnd(double mu) {
  return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
}

double truncgamma() {
  double c = MATH_PI_2;
  double X, gX;
  
  bool done = false;
  while(!done)
  {
    X = exprnd(1.0) * 2.0 + c;
    gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
    
    if(R::runif(0.0,1.0) <= gX) {
      done = true;
    }
  }
  
  return X;  
}

double randinvg(double mu) {
  // sampling
  double u = R::rnorm(0.0,1.0);
  double V = u*u;
  double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
  
  if(R::runif(0.0,1.0) > mu /(mu+out)) {    
    out = mu*mu / out; 
  }    
  return out;
}

double tinvgauss(double z, double t) {
  double X, u;
  double mu = 1.0/z;
  
  // Pick sampler
  if(mu > t) {
    // Sampler based on truncated gamma 
    // Algorithm 3 in the Windle (2013) PhD thesis, page 128
    while(1) {
      u = R::runif(0.0, 1.0);
      X = 1.0 / truncgamma();
      
      if ((double)std::log(u) < (-z*z*0.5*X)) {
        break;
      }
    }
  }  
  else {
    // Rejection sampler
    X = t + 1.0;
    while(X >= t) {
      X = randinvg(mu);
    }
  }    
  return X;
}

double samplepg(double z) {
  //  PG(b, z) = 0.25 * J*(b, z/2)
  z = (double)std::fabs((double)z) * 0.5;
  
  // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
  double t = MATH_2_PI;
  
  // Compute p, q and the ratio q / (q + p)
  // (derived from scratch; derivation is not in the original paper)
  double K = z*z/2.0 + MATH_PI2/8.0;
  double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
  double logK = (double)std::log(K);
  double Kt = K * t;
  double w = (double)std::sqrt(MATH_PI_2);
  
  double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
  double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
  double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
  double ratio = 1.0 / (1.0 + p_over_q); 
  
  double u, X;
  
  // Main sampling loop; page 130 of the Windle PhD thesis
  while(1) 
  {
    // Step 1: Sample X ? g(x|z)
    u = R::runif(0.0,1.0);
    if(u < ratio) {
      // truncated exponential
      X = t + exprnd(1.0)/K;
    }
    else {
      // truncated Inverse Gaussian
      X = tinvgauss(z, t);
    }
    
    // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
    int i = 1;
    double Sn = aterm(0, X, t);
    double U = R::runif(0.0,1.0) * Sn;
    int asgn = -1;
    bool even = false;
    
    while(1) 
    {
      Sn = Sn + asgn * aterm(i, X, t);
      
      // Accept if n is odd
      if(!even && (U <= Sn)) {
        X = X * 0.25;
        return X;
      }
      
      // Return to step 1 if n is even
      if(even && (U > Sn)) {
        break;
      }
      
      even = !even;
      asgn = -asgn;
      i++;
    }
  }
  return X;
}

double rpg(int n, double z){
  
  double x = 0;
  for(int i = 0; i < n; i++){
    x += samplepg(z);
  }
  
  return(x);
}

arma::vec samplePGvariables(arma::vec &Xbeta, arma::vec &n){
  
  int nsize = n.size();
  arma::vec Omega_vec(nsize);
  
  for(int i = 0; i < nsize; i++){
    
    Omega_vec[i] = rpg(n[i], Xbeta[i]);
    
  }
  
  return(Omega_vec);
}

///////////////  SAMPLE BETA PSI

arma::mat tXtOmegaX_betapsi(arma::mat &X, int Y, int X_centers, int ncov_psi, arma::vec Omega,
                            IntegerVector X_y_index, IntegerVector X_s_index){
  
  arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
  
  // year covariates times year covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
    
  }
  
  // spatial  covariates times spatial covariates
  for (int i = 0; (unsigned)i < X_s_index.size(); i++){
    
    XtOmegaX2(Y + X_s_index[i] - 1, Y + X_s_index[i] - 1) += Omega[i];
    
  }
  
  // year covariates times standard covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    for(int j = 0; j < ncov_psi; j++){
      
      XtOmegaX2(Y + X_centers + j, X_y_index[i] - 1) +=  X(i, Y + X_centers + j) * Omega[i];
      
    }
  }
  
  for (int i = 0; i < Y; i++) {
    for (int j = 0; j < ncov_psi; j++){
      XtOmegaX2(i, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i);
    }
  }
  // for (int i = 1; i <= Y; i++) {
  //   for (int j = 0; j < ncov_psi; j++){
  //     XtOmegaX2(i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i - 1);
  //   }
  // }
  
  // spatial  covariates times year covariates
  if(Y > 0){
    for (int i = 0; (unsigned)i < X_s_index.size(); i++){
      
      XtOmegaX2(X_y_index[i] - 1, Y + X_s_index[i] - 1) += Omega[i];
      
    }
  }
  
  for (int i = 0; i < Y; i++) {
    for (int j = 0; j < X_centers; j++){
      XtOmegaX2(Y + j, i) = XtOmegaX2(i, Y + j);
    }
  }
  // for (int i = 1; i <= Y; i++) {
  //   for (int j = 0; j < X_centers; j++){
  //     XtOmegaX2(Y + j, i - 1) = XtOmegaX2(i - 1, Y + j);
  //   }
  // }
  
  // spatial covariates times standard covariates
  for (int i = 0; (unsigned)i < X_s_index.size(); i++){
    
    for(int j = 0; j < ncov_psi; j++){
      
      XtOmegaX2(Y + X_centers + j, Y + X_s_index[i] - 1) +=  X(i, Y + X_centers + j) * Omega[i];
      
    }
  }
  
  for (int i = 1; i <= X_centers; i++) {
    
    for (int j = 0; j < ncov_psi; j++){
      XtOmegaX2(Y + i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, Y + i - 1);
    }
    
  }
  
  // standard covariates times standard covariates 
  
  for(int i = 0; i < ncov_psi; i++){
    for (int j = 0; j <= i; j++) {
      // arma::vec firstProduct = Omega % X.col(Y + X_centers + i);
      // arma::vec secondProduct = firstProduct % X.col(Y + X_centers + j);
      // XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = sum(secondProduct);
      for(int l = 0; l < Omega.size(); l++){
        XtOmegaX2(Y + X_centers + i, Y + X_centers + j) += Omega[l] * X(l, Y + X_centers + i) * X(l, Y + X_centers + j);
      }
    }
  }
  
  for (int i = 0; i < (ncov_psi- 1); i++) {
    for (int j = i; j < ncov_psi; j++) {
      XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = 
        XtOmegaX2(Y + X_centers + j, Y + X_centers + i);
    }
  }
  
  return(XtOmegaX2);
}

arma::vec XtransposeK_betapsi(arma::mat &X, IntegerVector X_y_index, 
                              IntegerVector X_s_index, arma::vec &k, 
                              int Y, int centers, int ncov){
  
  
  arma::vec Xk = arma::zeros(X.n_cols);
  
  for(int i = 0; i < X_y_index.size(); i++){
    Xk(X_y_index[i] - 1) += k[i];
  }
  for(int i = 0; i < X_s_index.size(); i++){
    Xk(Y + X_s_index[i] - 1) += k[i];
  }
  for(int i = 0; i < ncov; i++){
    Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
  }
  
  return(Xk);
}

arma::vec sampleBetaCoef_betapsi(arma::mat &X, arma::mat &invB, arma::vec &b, 
                                 arma::vec &k, arma::mat XtOmegaX,
                                 IntegerVector X_y_index, IntegerVector X_s_index,
                                 int Y, int centers, int ncov){
  
  // arma::mat tX = arma::trans(X);
  
  arma::mat tXk = XtransposeK_betapsi(X, X_y_index, X_s_index, k, Y, centers, ncov);
  
  arma::mat Lambda_B = XtOmegaX + invB;
  arma::vec mu_B = tXk + invB * b;
  
  arma::mat L = arma::trans(arma::chol(Lambda_B));
  arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
  
  // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
  
  arma::vec z = arma::randn(invB.n_cols);
  arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
  
  arma::vec result = v + alpha;
  
  return(result);
}

// [[Rcpp::export]]
List sampleBetaPsi(arma::vec beta,
                   arma::vec a_s,
                   arma::mat &X,
                   arma::vec Xbeta,
                   arma::vec b, 
                   arma::mat invB, 
                   arma::vec n, 
                   arma::vec k,
                   int Y, 
                   int X_centers,
                   int ncov_psi,
                   int numTimeSpaceCov,
                   IntegerVector X_y_index,
                   IntegerVector X_s_index){
  
  arma::vec Xbetaas = Xbeta + a_s;
  arma::vec Omega = samplePGvariables(Xbetaas, n);
  
  arma::mat XtOmegaX = tXtOmegaX_betapsi(X, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
                                         X_y_index, X_s_index);
  
  arma::vec knew = k - Omega % a_s;
  
  arma::vec betaout = sampleBetaCoef_betapsi(X, invB, b, knew, XtOmegaX,
                                             X_y_index, X_s_index, Y, X_centers, 
                                             numTimeSpaceCov + ncov_psi);
  
  return(List::create(_["beta"] = betaout,
                      _["Omega"] = Omega));
}

///////////////  SAMPLE BETA P

arma::mat XtOmegaX_betap(arma::mat &X, int p_intercepts, int ncov_p, arma::vec Omega,
                         IntegerVector X_y_index){
  
  arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
  
  
  // year covariates times year covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
    
  }
  
  // year covariates times standard covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    for(int j = 0; j < ncov_p; j++){
      
      XtOmegaX2(p_intercepts + j, X_y_index[i] - 1) +=  X(i, p_intercepts + j) * Omega[i];
      
    }
  }
  
  for (int i = 0; i < p_intercepts; i++) {
    for (int j = 0; j < ncov_p; j++){
      XtOmegaX2(i, p_intercepts + j) = XtOmegaX2(p_intercepts + j, i);
    }
  }
  
  // standard covariates times standard covariates 
  
  for(int i = 0; i < ncov_p; i++){
    for (int j = 0; j <= i; j++) {
      for(int l = 0; l < Omega.size(); l++){
        XtOmegaX2(p_intercepts + i, p_intercepts + j) += Omega[l] * 
          X(l, p_intercepts + i) * X(l, p_intercepts + j);
      }
    }
  }
  
  for (int i = 0; i < (ncov_p- 1); i++) {
    for (int j = i; j < ncov_p; j++) {
      XtOmegaX2(p_intercepts + i, p_intercepts + j) = 
        XtOmegaX2(p_intercepts + j, p_intercepts + i);
    }
  }
  
  
  
  return(XtOmegaX2);
}

arma::vec XtransposeK_betap(arma::mat &X, IntegerVector X_y_index, 
                            arma::vec &k, int p_intercepts, int ncov){
  
  
  arma::vec Xk = arma::zeros(X.n_cols);
  
  for(int i = 0; i < X_y_index.size(); i++){
    Xk(X_y_index[i] - 1) += k[i];
  }
  for(int i = 0; i < ncov; i++){
    Xk(p_intercepts + i) = as_scalar(k.t() * X.col(p_intercepts + i));
  }
  
  return(Xk);
}

arma::vec sampleBetaCoef_betap(arma::mat &X, arma::mat &invB, arma::vec &b, 
                               arma::vec &k, arma::mat XtOmegaX,
                               IntegerVector X_y_index,
                               int p_intercepts, int ncov){
  
  // arma::mat tX = arma::trans(X);
  
  arma::mat tXk = XtransposeK_betap(X, X_y_index, k, p_intercepts, ncov);
  
  arma::mat Lambda_B = XtOmegaX + invB;
  arma::vec mu_B = tXk + invB * b;
  
  arma::mat L = arma::trans(arma::chol(Lambda_B));
  arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
  arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
  
  // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
  
  arma::vec z = arma::randn(invB.n_cols);
  arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
  
  arma::vec result = v + alpha;
  
  return(result);
}

// [[Rcpp::export]]
arma::vec sampleBetaP(arma::vec beta,
                      arma::mat &X, 
                      arma::vec Xbeta,
                      arma::vec b, 
                      arma::mat invB,
                      int ncov_p,
                      int p_intercepts,
                      IntegerVector X_y_index,
                      arma::vec n, 
                      arma::vec k){
  
  arma::vec Omega = samplePGvariables(Xbeta, n);
  
  arma::mat XtOmegaX = XtOmegaX_betap(X, p_intercepts, ncov_p, Omega,
                                      X_y_index);
  
  beta = sampleBetaCoef_betap(X, invB, b, k, XtOmegaX, X_y_index, p_intercepts, ncov_p);
  
  return(beta);
}

// GAUSSIAN PROCESS FUNCTIONS

double k_cpp(double x1, double x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-(x1-x2)*(x1-x2)/(2*pow(l,2)));
  // return 1;
}

// [[Rcpp::export]]
arma::mat K(arma::vec x1, arma::vec x2, double a, double l){
  arma::mat res(x1.size(), x2.size());
  
  for(int i = 0; (unsigned)i < x1.size(); i++){
    for(int j = 0; (unsigned)j < x2.size(); j++){
      res(i,j) = k_cpp(x1[i],x2[j], a, l);
    }  
  }
  
  return res;
}

double k2_cpp(arma::rowvec x1, arma::rowvec x2, double a, double l){
  // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
  return a*exp(-( pow(x1[0]-x2[0], 2) + pow(x1[1]-x2[1], 2) ) /(2*pow(l,2)));
}

// [[Rcpp::export]]
arma::mat K2(arma::mat x1, arma::mat x2, double a, double l){
  arma::mat res(x1.n_rows, x2.n_rows);
  
  for(int i = 0; (unsigned)i < x1.n_rows; i++){
    for(int j = 0; (unsigned)j < x2.n_rows; j++){
      res(i,j) = k2_cpp(x1.row(i),x2.row(j), a, l);
    }  
  }
  
  return res;
}

// SAMPLE SITE-SPECIFIC RANDOM EFFECTS

// [[Rcpp::export]]
arma::vec sample_eps_cpp(arma::vec k_s, arma::vec sites,
                         arma::vec &c_psi,
                         arma::vec k,
                         arma::vec z, arma::vec Omega,
                         double sigma_a){
  
  arma::vec b_psi = k - Omega % c_psi;
  
  arma::vec a_s = arma::zeros(k_s.size());
  
  int index_site = 0;
  
  for(int i = 0; (unsigned)i < sites.size(); i++){
    
    int site = sites[i];
    
    int l = 1;
    
    // find rows associated with current site
    if((unsigned)i != (sites.size() - 1)){
      
      while(k_s[index_site + l] == site){
        l += 1;
      }
      
    } else {
      
      l = k_s.size() - index_site;
      
    }
    
    IntegerVector indexes_site(l);
    for(int j = 0; j < l; j++){
      indexes_site[j] = index_site + j;
    }
    index_site += l;
    
    // arma::mat data_psi_site(l, X_psi.n_cols);
    // arma::vec z_site = arma::zeros(l);
    arma::vec w_site = arma::zeros(l);
    // arma::vec c_site = arma::zeros(l);
    arma::vec b_site = arma::zeros(l);
    for(int j = 0; j < l; j++){
      // data_psi_site.row(j) = X_psi.row(indexes_site[j]);
      // z_site[j] = z[indexes_site[j]];
      w_site[j] = Omega[indexes_site[j]];
      // c_site[j] = c_psi[indexes_site[j]];
      b_site[j] = b_psi[indexes_site[j]];
    }
    
    // sample x
    // arma::vec k_site = (z_site - .5);
    // arma::vec n_site = arma::ones(l);
    // arma::vec c_site = data_psi_site * beta_psi;
    
    double a = 1 / (sum(w_site) + 1 / (sigma_a*sigma_a));
    // double b = sum(k_site - w_site % c_site);
    double b = sum(b_site);
    
    double x = R::rnorm(b * a, sqrt(a));
    
    for(int j = 0; j < l; j++){
      a_s[indexes_site[j]] = x;
    }
    
  }
  
  return(a_s);
}

// CREATE GRID

// [[Rcpp::export]]
bool isPointInBandRight(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
  
  for(int k = 0; k < X_tilde.n_rows; k++){
    
    if((X_tilde(k,1) < y_grid[j + 1]) & (X_tilde(k,1) > y_grid[j - 1])){
      if(X_tilde(k,0) < x_grid[i + 1]){
        return(true);
      }
    } 
    
  }
  
  return(false);
}

// [[Rcpp::export]]
bool isPointInBandLeft(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j) {
  
  for(int k = 0; k < X_tilde.n_rows; k++){
    
    if((X_tilde(k,1) < y_grid[j + 1]) & (X_tilde(k,1) > y_grid[j - 1])){
      if(X_tilde(k,0) > x_grid[i - 1]){
        return(true);
      }
    } 
    
  }
  
  return(false);
}

// [[Rcpp::export]]
bool isPointInBandUp(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
  
  for(int k = 0; k < X_tilde.n_rows; k++){
    
    if((X_tilde(k,0) < x_grid[i + 1]) & (X_tilde(k,0) > x_grid[i - 1])){
      if(X_tilde(k,1) > y_grid[j-1]){
        return(true);
      }
    }
    
  }
  
  return(false);
  
}

// [[Rcpp::export]]
bool isPointInBandDown(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
  
  for(int k = 0; k < X_tilde.n_rows; k++){
    
    if((X_tilde(k,0) < x_grid[i + 1]) & (X_tilde(k,0) > x_grid[i - 1])){
      if(X_tilde(k,1) < y_grid[j+1]){
        return(true);
      }
    }
    
  }
  
  return(false);
  
}

// [[Rcpp::export]]
IntegerVector findClosestPoint(arma::mat XY_sp, arma::mat X_tilde){
  
  IntegerVector closestPoint(XY_sp.n_rows);
  
  for(int k = 0; k < XY_sp.n_rows; k++){
    
    double newDistance = 0;
    double minDistance = exp(10);
    int bestIndex = 0;
    
    for(int i = 0; i < X_tilde.n_rows; i++){
      newDistance = pow(X_tilde(i, 0) - XY_sp(k, 0), 2) + pow(X_tilde(i, 1) - XY_sp(k, 1), 2);
      
      if(newDistance < minDistance){
        minDistance = newDistance;
        bestIndex = i + 1;
      }
    }
    
    closestPoint[k] = bestIndex;
    
  }
  
  return(closestPoint);
}

//// HYPERPARAMETER SAMPLING FUNCTIONS

// [[Rcpp::export]]
arma::mat matrixProductXtOmegaX_year(int Y, arma::vec Omega,
                                     arma::vec X_y_index){
  
  arma::mat XtOmegaX2 = arma::zeros(Y, Y);
  
  // year covariates times year covariates
  for (int i = 0; (unsigned)i < X_y_index.size(); i++){
    
    XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
    
  }
  
  return(XtOmegaX2);
}

// [[Rcpp::export]]
arma::mat matrixProductXtOmegaX_spatial(int X_centers, arma::vec Omega,
                                        arma::vec X_s_index){
  
  arma::mat XtOmegaX2 = arma::zeros(X_centers, X_centers);
  
  // spatial  covariates times spatial covariates
  for (int i = 0; (unsigned)i < X_s_index.size(); i++){
    
    XtOmegaX2(X_s_index[i] - 1, X_s_index[i] - 1) += Omega[i];
    
  }
  
  return(XtOmegaX2);
}

// [[Rcpp::export]]
arma::vec XpsiYz(arma::vec X_y_index, arma::vec z, int Y){
  
  arma::vec out = arma::zeros(Y);
  
  for(int i = 0; i < X_y_index.size(); i++){
    out(X_y_index[i] - 1) += z[i];
  }
  
  return(out);
}

// [[Rcpp::export]]
arma::vec XpsinoYbetaz(arma::vec X_s_index, int X_centers, int ncov,
                       arma::mat &X_cov, arma::vec& beta){
  
  // X_psi[,-(1:Y),drop=F] %*% beta_psi[-(1:Y)]
  
  return(X_cov * beta);
  // arma::vec out = arma::zeros(X.n_cols);
  
  // for(int i = 0; i < X_s_index.size(); i++){
  //   Xk(Y + X_s_index[i] - 1) += k[i];
  // }
  // for(int i = 0; i < ncov; i++){
  //   Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
  // }
  
}
