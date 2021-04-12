#include <RcppArmadilloExtensions/sample.h>

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

Rcpp::NumericVector arma_setdiff(arma::uvec& x, arma::uvec& y){
  
  x = arma::unique(x);
  y = arma::unique(y);
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uvec q1 = arma::find(x == y[j]);
    x.shed_row(q1(0));
  }
  
  Rcpp::NumericVector x2 = Rcpp::wrap(x);
  x2.attr("dim") = R_NilValue;
  return x2;
}

IntegerVector stl_sort(IntegerVector x) {
  IntegerVector y = clone(x);
  std::sort(y.begin(), y.end());
  return y;
}

IntegerVector setdiffna(IntegerVector S_c){
  IntegerVector S_cnew = stl_sort(unique(S_c));
  IntegerVector NAvec = IntegerVector::create(NA_INTEGER);
  arma::uvec temp1 = as<arma::uvec>(S_cnew);
  arma::uvec temp2 = as<arma::uvec>(NAvec);
  IntegerVector temp = as<IntegerVector>(arma_setdiff(temp1, temp2));
  return temp;
}

///////// SAMPLE LATENT OCCUPANCIES

//' Sample the occupancies
//' 
//' @param psi Probabilities of presence.
//' @param p Probabilities of detection.
//' @param k_s Correspondance.
//' @export
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

// ///////// GOODNESS OF FIT
// 
// // [[Rcpp::export]]
// arma::vec simulateDetections(arma::vec p, arma::vec z_all){
//   
//   arma::vec y = arma::zeros(p.size());
//   
//   // this loops over p
//   for(int i = 0; (unsigned)i < y.size(); i++){
//     
//     if(z_all(i) == 1){
//       
//       y[i] = R::rbinom(1, p[i]);
//       
//     } 
//     
//   }
//   
//   return(y);
// }  
// 
// arma::vec logit(arma::vec x){
//   return(1 / (1 + exp(-x)));  
// }
// 
// // [[Rcpp::export]]
// arma::vec computeYearEffect(int Y, arma::vec a_s_unique, arma::vec beta_psi){
//   
//   arma::vec yearEffect = arma::zeros(Y);
//   
//   for(int y = 0; (unsigned)y < Y; y++){
//     
//     arma::vec occProbs = logit(beta_psi[y] + a_s_unique);
//     
//     yearEffect[y] = mean(occProbs);
//   }
//   
//   return(yearEffect);
// }  
// 
// //////// SAMPLE POLYA GAMMA VARIABLES
// 
// double aterm(int n, double x, double t) {
//   double f = 0;
//   if(x <= t) {
//     f = MATH_LOG_PI + (double)std::log(n + 0.5) + 1.5*(MATH_LOG_2_PI- (double)std::log(x)) - 2*(n + 0.5)*(n + 0.5)/x;
//   }
//   else {
//     f = MATH_LOG_PI + (double)std::log(n + 0.5) - x * MATH_PI2_2 * (n + 0.5)*(n + 0.5);
//   }    
//   return (double)exp(f);
// }
// 
// double exprnd(double mu) {
//   return -mu * (double)std::log(1.0 - (double)R::runif(0.0,1.0));
// }
// 
// double truncgamma() {
//   double c = MATH_PI_2;
//   double X, gX;
//   
//   bool done = false;
//   while(!done)
//   {
//     X = exprnd(1.0) * 2.0 + c;
//     gX = MATH_SQRT_PI_2 / (double)std::sqrt(X);
//     
//     if(R::runif(0.0,1.0) <= gX) {
//       done = true;
//     }
//   }
//   
//   return X;  
// }
// 
// double randinvg(double mu) {
//   // sampling
//   double u = R::rnorm(0.0,1.0);
//   double V = u*u;
//   double out = mu + 0.5*mu * ( mu*V - (double)std::sqrt(4.0*mu*V + mu*mu * V*V) );
//   
//   if(R::runif(0.0,1.0) > mu /(mu+out)) {    
//     out = mu*mu / out; 
//   }    
//   return out;
// }
// 
// double tinvgauss(double z, double t) {
//   double X, u;
//   double mu = 1.0/z;
//   
//   // Pick sampler
//   if(mu > t) {
//     // Sampler based on truncated gamma 
//     // Algorithm 3 in the Windle (2013) PhD thesis, page 128
//     while(1) {
//       u = R::runif(0.0, 1.0);
//       X = 1.0 / truncgamma();
//       
//       if ((double)std::log(u) < (-z*z*0.5*X)) {
//         break;
//       }
//     }
//   }  
//   else {
//     // Rejection sampler
//     X = t + 1.0;
//     while(X >= t) {
//       X = randinvg(mu);
//     }
//   }    
//   return X;
// }
// 
// double samplepg(double z) {
//   //  PG(b, z) = 0.25 * J*(b, z/2)
//   z = (double)std::fabs((double)z) * 0.5;
//   
//   // Point on the intersection IL = [0, 4/ log 3] and IR = [(log 3)/pi^2, \infty)
//   double t = MATH_2_PI;
//   
//   // Compute p, q and the ratio q / (q + p)
//   // (derived from scratch; derivation is not in the original paper)
//   double K = z*z/2.0 + MATH_PI2/8.0;
//   double logA = (double)std::log(4.0) - MATH_LOG_PI - z;
//   double logK = (double)std::log(K);
//   double Kt = K * t;
//   double w = (double)std::sqrt(MATH_PI_2);
//   
//   double logf1 = logA + R::pnorm(w*(t*z - 1),0.0,1.0,1,1) + logK + Kt;
//   double logf2 = logA + 2*z + R::pnorm(-w*(t*z+1),0.0,1.0,1,1) + logK + Kt;
//   double p_over_q = (double)std::exp(logf1) + (double)std::exp(logf2);
//   double ratio = 1.0 / (1.0 + p_over_q); 
//   
//   double u, X;
//   
//   // Main sampling loop; page 130 of the Windle PhD thesis
//   while(1) 
//   {
//     // Step 1: Sample X ? g(x|z)
//     u = R::runif(0.0,1.0);
//     if(u < ratio) {
//       // truncated exponential
//       X = t + exprnd(1.0)/K;
//     }
//     else {
//       // truncated Inverse Gaussian
//       X = tinvgauss(z, t);
//     }
//     
//     // Step 2: Iteratively calculate Sn(X|z), starting at S1(X|z), until U ? Sn(X|z) for an odd n or U > Sn(X|z) for an even n
//     int i = 1;
//     double Sn = aterm(0, X, t);
//     double U = R::runif(0.0,1.0) * Sn;
//     int asgn = -1;
//     bool even = false;
//     
//     while(1) 
//     {
//       Sn = Sn + asgn * aterm(i, X, t);
//       
//       // Accept if n is odd
//       if(!even && (U <= Sn)) {
//         X = X * 0.25;
//         return X;
//       }
//       
//       // Return to step 1 if n is even
//       if(even && (U > Sn)) {
//         break;
//       }
//       
//       even = !even;
//       asgn = -asgn;
//       i++;
//     }
//   }
//   return X;
// }
// 
// // [[Rcpp::export]]
// double rpg(int n, double z){
//   
//   double x = 0;
//   for(int i = 0; i < n; i++){
//     x += samplepg(z);
//   }
//   
//   return(x);
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_Omega_cpp(arma::mat &X, arma::vec &beta, arma::vec &n){
//   
//   int nsize = n.size();
//   arma::vec Omega_vec(nsize);
//   
//   arma::vec Xbeta = X * beta;
//   
//   for(int i = 0; i < nsize; i++){
//     
//     Omega_vec[i] = rpg(n[i], Xbeta[i]);
//     
//   }
//   
//   return(Omega_vec);
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_Omega_cpp_noXb(arma::vec &Xbeta, arma::vec &n){
//   
//   int nsize = n.size();
//   arma::vec Omega_vec(nsize);
//   
//   for(int i = 0; i < nsize; i++){
//     
//     Omega_vec[i] = rpg(n[i], Xbeta[i]);
//     
//   }
//   
//   return(Omega_vec);
// }
// 
// struct Sampleomega : public Worker {
//   
//   // inputs
//   RMatrix<double> X;
//   RVector<double> beta;
//   RVector<double> n_trials;
//   
//   // output matrix to write to
//   RVector<double> Omega;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   Sampleomega(NumericMatrix X, NumericVector beta, NumericVector n_trials, NumericVector Omega)
//     : X(X), beta(beta), n_trials(n_trials), Omega(Omega) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       RMatrix<double>::Row row1 = X.row(i);
//       // double myRes = row1 * beta;
//       double Xibeta = std::inner_product(row1.begin(),
//                                          row1.end(),
//                                          beta.begin(),
//                                          0.0);
//       
//       Omega[i] = rpg(n_trials[i], Xibeta);
//       
//     }
//   }
// }; 
// 
// // [[Rcpp::export]]
// NumericVector sampleOmegaParallel(NumericMatrix &X, NumericVector beta, NumericVector n_trials) {
//   
//   // allocate the matrix we will return
//   NumericVector Omega(X.nrow());
//   
//   // create the worker
//   Sampleomega sampleomega(X, beta, n_trials, Omega);
//   
//   // call it with parallelFor
//   parallelFor(0, X.nrow(), sampleomega);
//   
//   return Omega;
// }
// 
// struct Sampleomegaas : public Worker {
//   
//   // inputs
//   RMatrix<double> X;
//   RVector<double> beta;
//   RVector<double> as;
//   RVector<double> n_trials;
//   
//   // output matrix to write to
//   RVector<double> Omega;
//   
//   // initialize from Rcpp input and output matrixes (the RMatrix class
//   // can be automatically converted to from the Rcpp matrix type)
//   Sampleomegaas(NumericMatrix X, NumericVector beta, NumericVector as,
//                 NumericVector n_trials, NumericVector Omega)
//     : X(X), beta(beta), as(as), n_trials(n_trials), Omega(Omega) {}
//   
//   // function call operator that work for the specified range (begin/end)
//   void operator()(std::size_t begin, std::size_t end) {
//     for (std::size_t i = begin; i < end; i++) {
//       
//       RMatrix<double>::Row row1 = X.row(i);
//       // double myRes = row1 * beta;
//       double Xibeta = std::inner_product(row1.begin(),
//                                          row1.end(),
//                                          beta.begin(),
//                                          0.0);
//       Xibeta += as[i];
//       Omega[i] = rpg(n_trials[i], Xibeta);
//       
//     }
//   }
// };
// 
// // [[Rcpp::export]]
// NumericVector sampleOmegaParallelas(NumericMatrix &X, NumericVector beta, 
//                                     NumericVector as, NumericVector n_trials) {
//   
//   // allocate the matrix we will return
//   NumericVector Omega(X.nrow());
//   
//   // create the worker
//   Sampleomegaas sampleomega(X, beta, as, n_trials, Omega);
//   
//   // call it with parallelFor
//   parallelFor(0, X.nrow(), sampleomega);
//   
//   return Omega;
// }
// 
// ///////////// DISTRIBUTIONS
// 
// arma::vec mvrnormArma(arma::vec mu, arma::mat sigma) {
//   int ncols = sigma.n_cols;
//   arma::vec Y = arma::randn(ncols);
//   return mu + arma::chol(sigma) * Y;
// }
// 
// arma::vec mvrnormArmaQuick(arma::vec mu, arma::mat cholsigma) {
//   int ncols = cholsigma.n_cols;
//   arma::vec Y = arma::randn(ncols);
//   return mu + cholsigma * Y;
// }
// 
// // [[Rcpp::export]]
// double dt_cpp(double x, double nu, double mu, double sigma, bool returnLog){
//   
//   double ratio = exp(R::lgammafn(( nu + 1 ) / 2) - R::lgammafn( nu / 2 ));
//   double product = (x - mu) * (x - mu) / (sigma * sigma);
//   double num = pow(1 + (1 / nu) * product, - ( nu + 1 ) / 2);
//   double den = sigma * sqrt(M_PI * nu);
//   if(returnLog){
//     return log(ratio * num / den);
//   } else {
//     return ratio * num / den;
//   }
// }
// 
// ////////////// MATRIX PRODUCTS
// 
// // [[Rcpp::export]]
// arma::mat diagMatrixProd(arma::mat &X, arma::vec &D){
//   // this is slow
//   arma::mat result(X.n_rows, D.size());
//   for(int j = 0; (unsigned)j < result.n_cols; j++){
//     result.col(j) = X.col(j) * D(j);
//   }
//   
//   // RMatrix<double>::Row rowi = A.row(i);
//   
//   // out(i,1) = std::inner_product(rowi.begin(), rowi.end(), x.begin(), 0.0);
//   
//   return(result);
// }
// 
// // [[Rcpp::export]]
// arma::mat matrixProductXtOmegaX(arma::mat &X, int Y, int X_centers, int ncov_psi, arma::vec Omega,
//                                 IntegerVector X_y_index, IntegerVector X_s_index){
//   
//   arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
//   
//   // year covariates times year covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
//     
//   }
//   
//   // spatial  covariates times spatial covariates
//   for (int i = 0; (unsigned)i < X_s_index.size(); i++){
//     
//     XtOmegaX2(Y + X_s_index[i] - 1, Y + X_s_index[i] - 1) += Omega[i];
//     
//   }
//   
//   // year covariates times standard covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       XtOmegaX2(Y + X_centers + j, X_y_index[i] - 1) +=  X(i, Y + X_centers + j) * Omega[i];
//       
//     }
//   }
//   
//   for (int i = 0; i < Y; i++) {
//     for (int j = 0; j < ncov_psi; j++){
//       XtOmegaX2(i, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i);
//     }
//   }
//   // for (int i = 1; i <= Y; i++) {
//   //   for (int j = 0; j < ncov_psi; j++){
//   //     XtOmegaX2(i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i - 1);
//   //   }
//   // }
//   
//   // spatial  covariates times year covariates
//   if(Y > 0){
//     for (int i = 0; (unsigned)i < X_s_index.size(); i++){
//       
//       XtOmegaX2(X_y_index[i] - 1, Y + X_s_index[i] - 1) += Omega[i];
//       
//     }
//   }
//   
//   for (int i = 0; i < Y; i++) {
//     for (int j = 0; j < X_centers; j++){
//       XtOmegaX2(Y + j, i) = XtOmegaX2(i, Y + j);
//     }
//   }
//   // for (int i = 1; i <= Y; i++) {
//   //   for (int j = 0; j < X_centers; j++){
//   //     XtOmegaX2(Y + j, i - 1) = XtOmegaX2(i - 1, Y + j);
//   //   }
//   // }
//   
//   // spatial covariates times standard covariates
//   for (int i = 0; (unsigned)i < X_s_index.size(); i++){
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       XtOmegaX2(Y + X_centers + j, Y + X_s_index[i] - 1) +=  X(i, Y + X_centers + j) * Omega[i];
//       
//     }
//   }
//   
//   for (int i = 1; i <= X_centers; i++) {
//     
//     for (int j = 0; j < ncov_psi; j++){
//       XtOmegaX2(Y + i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, Y + i - 1);
//     }
//     
//   }
//   
//   // standard covariates times standard covariates 
//   
//   for(int i = 0; i < ncov_psi; i++){
//     for (int j = 0; j <= i; j++) {
//       // arma::vec firstProduct = Omega % X.col(Y + X_centers + i);
//       // arma::vec secondProduct = firstProduct % X.col(Y + X_centers + j);
//       // XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = sum(secondProduct);
//       for(int l = 0; l < Omega.size(); l++){
//         XtOmegaX2(Y + X_centers + i, Y + X_centers + j) += Omega[l] * X(l, Y + X_centers + i) * X(l, Y + X_centers + j);
//       }
//     }
//   }
//   
//   for (int i = 0; i < (ncov_psi- 1); i++) {
//     for (int j = i; j < ncov_psi; j++) {
//       XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = 
//         XtOmegaX2(Y + X_centers + j, Y + X_centers + i);
//     }
//   }
//   
//   return(XtOmegaX2);
// }
// 
// // [[Rcpp::export]]
// arma::mat matrixProductXtOmegaX_SoR(arma::mat &X, int Y, int X_centers, int ncov_psi, arma::vec Omega,
//                                     IntegerVector X_y_index, arma::mat X_s_index){
//   
//   arma::mat XtOmegaX2 = arma::zeros(X.n_cols, X.n_cols);
//   
//   int maxPoints = X_s_index.n_cols;
//   
//   // year covariates times year covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
//     
//   }
//   
//   // spatial  covariates times spatial covariates
//   for(int l = 0; l < maxPoints; l++){
//     
//     for(int l2 = 0; l2 < maxPoints; l2++){
//       
//       for (int i = 0; (unsigned)i < Omega.size(); i++){
//         
//         XtOmegaX2(Y + X_s_index(i,l) - 1, Y + X_s_index(i,l2) - 1) += X(i, Y + X_s_index(i, l) - 1) * 
//           Omega[i] * X(i, Y + X_s_index(i, l2) - 1);    
//         
//       }
//       
//       // XtOmegaX2(Y + X_s_index[i,l2] - 1, Y + X_s_index[i,l] - 1) = 
//       // XtOmegaX2(Y + X_s_index[i,l] - 1, Y + X_s_index[i,l2] - 1);
//       
//     }
//     
//   }
//   
//   // year covariates times standard covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       XtOmegaX2(Y + X_centers + j, X_y_index[i] - 1) +=  X(i, Y + X_centers + j) * Omega[i];
//       
//     }
//   }
//   
//   for (int i = 1; i <= Y; i++) {
//     
//     for (int j = 0; j < ncov_psi; j++){
//       XtOmegaX2(i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, i - 1);
//     }
//     
//   }
//   
//   // spatial  covariates times year covariates
//   for (int i = 0; (unsigned)i < Omega.size(); i++){
//     
//     for(int l = 0; l < maxPoints; l++){
//       
//       XtOmegaX2(X_y_index[i] - 1, Y + X_s_index(i,l) - 1) += X(i, Y + X_s_index(i,l) - 1) * Omega[i];
//       
//     }
//     
//   }
//   
//   for (int i = 1; i <= Y; i++) {
//     
//     for (int j = 0; j < X_centers; j++){
//       XtOmegaX2(Y + j, i - 1) = XtOmegaX2(i - 1, Y + j);
//     }
//     
//   }
//   
//   // spatial covariates times standard covariates
//   for (int i = 0; (unsigned)i < Omega.size(); i++){
//     
//     for(int j = 0; j < ncov_psi; j++){
//       
//       for(int l = 0; l < maxPoints; l++){
//         
//         XtOmegaX2(Y + X_centers + j, Y + X_s_index(i,l) - 1) +=  
//           X(i, Y + X_centers + j) * X(i, Y + X_s_index(i,l) - 1) * Omega[i];
//         
//       }
//     }
//   }
//   
//   for (int i = 1; i <= X_centers; i++) {
//     
//     for (int j = 0; j < ncov_psi; j++){
//       XtOmegaX2(Y + i - 1, Y + X_centers + j) = XtOmegaX2(Y + X_centers + j, Y + i - 1);
//     }
//     
//   }
//   
//   // standard covariates times standard covariates 
//   
//   for(int i = 0; i < ncov_psi; i++){
//     for (int j = 0; j <= i; j++) {
//       // arma::vec firstProduct = Omega % X.col(Y + X_centers + i);
//       // arma::vec secondProduct = firstProduct % X.col(Y + X_centers + j);
//       // XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = sum(secondProduct);
//       for(int l = 0; l < Omega.size(); l++){
//         XtOmegaX2(Y + X_centers + i, Y + X_centers + j) += Omega[l] * X(l, Y + X_centers + i) * X(l, Y + X_centers + j);
//       }
//     }
//   }
//   
//   for (int i = 0; i < (ncov_psi- 1); i++) {
//     for (int j = i; j < ncov_psi; j++) {
//       XtOmegaX2(Y + X_centers + i, Y + X_centers + j) = 
//         XtOmegaX2(Y + X_centers + j, Y + X_centers + i);
//     }
//   }
//   
//   return(XtOmegaX2);
// }
// 
// ///////////////  SAMPLE COEFFICIENTS
// 
// // [[Rcpp::export]]
// arma::vec sample_beta_cpp_fast(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                arma::vec &k, arma::mat XtOmegaX){
//   
//   // arma::mat tX = ;
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = arma::trans(X) * k + invB * b;
//   
//   // arma::mat L = arma::trans(arma::chol(XtOmegaX + invB));
//   // arma::vec tmp = arma::solve(arma::trimatl(L), tX * k + invB * b);
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   return(result);
// }
// 
// // [[Rcpp::export]]
// arma::vec XtransposeKarma(arma::mat &X, IntegerVector X_y_index, 
//                           IntegerVector X_s_index, arma::vec &k, 
//                           int Y, int centers, int ncov){
//   
//   
//   arma::vec Xk = arma::zeros(X.n_cols);
//   
//   for(int i = 0; i < X_y_index.size(); i++){
//     Xk(X_y_index[i] - 1) += k[i];
//   }
//   for(int i = 0; i < X_s_index.size(); i++){
//     Xk(Y + X_s_index[i] - 1) += k[i];
//   }
//   for(int i = 0; i < ncov; i++){
//     Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
//   }
//   
//   return(Xk);
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_beta_cpp_fast_sparse(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                       arma::vec &k, arma::mat XtOmegaX,
//                                       IntegerVector X_y_index, IntegerVector X_s_index,
//                                       int Y, int centers, int ncov){
//   
//   // arma::mat tX = arma::trans(X);
//   
//   arma::mat tXk = XtransposeKarma(X, X_y_index, X_s_index, k, Y, centers, ncov);
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = tXk + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   return(result);
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_beta_cpp_fast_sparse_constrained(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                                   arma::vec &k, arma::mat XtOmegaX,
//                                                   IntegerVector X_y_index, IntegerVector X_s_index,
//                                                   int Y, int centers, int ncov){
//   
//   // arma::mat tX = arma::trans(X);
//   
//   arma::mat tXk = XtransposeKarma(X, X_y_index, X_s_index, k, Y, centers, ncov);
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = tXk + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   // constraint
//   
//   // arma::mat invQ = Lambda_B;
//   arma::mat invQ = arma::inv(Lambda_B);
//   
//   // sQ <- sum(Lambda_B[Y + 1:X_centers, Y + 1:X_centers])
//   double sQ = 0;
//   for(int l = 0; l < centers; l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       sQ += invQ(Y + l, Y + l2);
//     }
//   }
//   
//   // columns <- apply(Lambda_B[,Y + 1:X_centers], 1, sum) / sQ
//   arma::vec columns = arma::zeros(b.size());
//   for(int l = 0; l < b.size(); l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       columns[l] += invQ(l, Y + l2);
//     }
//     columns[l] = columns[l] / sQ;
//   }
//   
//   // allTerm <- sapply(1:length(x), function(i){
//   //sum(columns[i] * x[Y + 1:X_centers])
//   //})
//   arma::vec secondTerm = arma::zeros(b.size());
//   for(int l = 0; l < b.size(); l++){
//     for(int l2 = 0; l2 < centers; l2++){
//       secondTerm[l] += columns[l] * result[Y + l2];
//     }
//   }
//   
//   // Rcout << columns << " - " << sQ << std::endl;
//   
//   result = result - secondTerm;
//   
//   return(result);
// }
// 
// // [[Rcpp::export]]
// arma::vec XtransposeKarma_SoR(arma::mat &X, IntegerVector X_y_index, 
//                               arma::mat &X_s_index, arma::vec &k, 
//                               int Y, int centers, int ncov){
//   
//   
//   arma::vec Xk = arma::zeros(X.n_cols);
//   
//   for(int i = 0; i < X_y_index.size(); i++){
//     Xk(X_y_index[i] - 1) += k[i];
//   }
//   
//   for(int i = 0; i < k.size(); i++){
//     for(int l = 0; l < X_s_index.n_cols; l++){
//       Xk(Y + X_s_index(i,l) - 1) += X(i, Y + X_s_index(i,l) - 1) * k[i];
//     }
//   }
//   
//   for(int i = 0; i < ncov; i++){
//     Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
//   }
//   
//   return(Xk);
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_beta_cpp_fast_sparse_SoR(arma::mat &X, arma::mat &invB, arma::vec &b, 
//                                           arma::vec &k, arma::mat XtOmegaX,
//                                           IntegerVector X_y_index, arma::mat &X_s_index,
//                                           int Y, int centers, int ncov){
//   
//   // arma::mat tX = arma::trans(X);
//   
//   arma::mat tXk = XtransposeKarma_SoR(X, X_y_index, X_s_index, k, Y, centers, ncov);
//   
//   arma::mat Lambda_B = XtOmegaX + invB;
//   arma::vec mu_B = tXk + invB * b;
//   
//   arma::mat L = arma::trans(arma::chol(Lambda_B));
//   arma::vec tmp = arma::solve(arma::trimatl(L), mu_B);
//   arma::vec alpha = arma::solve(arma::trimatu(arma::trans(L)),tmp);
//   
//   // arma::vec result = mvrnormArmaQuick(alpha, arma::trans(arma::inv(arma::trimatl(L))));
//   
//   arma::vec z = arma::randn(invB.n_cols);
//   arma::vec v = arma::solve(arma::trimatu(arma::trans(L)), z);
//   
//   arma::vec result = v + alpha;
//   
//   return(result);
// }
// 
// 
// // [[Rcpp::export]]
// List sampler_beta(arma::vec beta,
//                   arma::vec a_s,
//                   arma::mat &X, 
//                   arma::vec b, 
//                   arma::mat invB, 
//                   arma::vec n, 
//                   arma::vec k,
//                   int Y, 
//                   int X_centers,
//                   int ncov_psi,
//                   int numTimeSpaceCov,
//                   IntegerVector X_y_index,
//                   IntegerVector X_s_index){
//   
//   arma::vec Xbeta = X * beta + a_s;
//   arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX(X, Y, X_centers, ncov_psi + numTimeSpaceCov, Omega,
//                                              X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % a_s;
//   
//   beta = sample_beta_cpp_fast_sparse(X, invB, b, knew, XtOmegaX,
//                                      X_y_index, X_s_index, Y, X_centers, ncov_psi + numTimeSpaceCov);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }
// 
// // [[Rcpp::export]]
// List sampler_beta_parallel(SEXP beta,
//                            SEXP a_s,
//                            SEXP &X, 
//                            arma::vec b, 
//                            arma::mat invB, 
//                            SEXP n, 
//                            arma::vec k,
//                            int Y, 
//                            int X_centers,
//                            int ncov_psi,
//                            int numTimeSpaceCov,
//                            IntegerVector X_y_index,
//                            IntegerVector X_s_index){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   NumericVector as_num = as<NumericVector>(a_s);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallelas(X_num, beta_num, as_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   arma::vec as_arma(as_num.begin(), as_num.size(), false);
//   
//   // arma::vec Xbeta = X * beta + a_s;
//   // arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX(X_arma, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
//                                              X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % as_arma;
//   
//   // arma::vec betaout = sample_beta_cpp_fast_sparse(X_arma, invB, b, knew, XtOmegaX,
//   // X_y_index, X_s_index, Y, X_centers, numTimeSpaceCov + ncov_psi);
//   arma::vec betaout = sample_beta_cpp_fast_sparse(X_arma, invB, b, knew, XtOmegaX,
//                                                   X_y_index, X_s_index, Y, X_centers, 
//                                                   numTimeSpaceCov + ncov_psi);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }
// 
// // [[Rcpp::export]]
// List sampler_beta_parallel_SoR(SEXP beta,
//                                SEXP a_s,
//                                SEXP &X, 
//                                arma::vec b, 
//                                arma::mat invB, 
//                                SEXP n, 
//                                arma::vec k,
//                                int Y, 
//                                int X_centers,
//                                int ncov_psi,
//                                int numTimeSpaceCov,
//                                IntegerVector X_y_index,
//                                arma::mat &X_s_index){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   NumericVector as_num = as<NumericVector>(a_s);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallelas(X_num, beta_num, as_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   arma::vec as_arma(as_num.begin(), as_num.size(), false);
//   
//   // arma::vec Xbeta = X * beta + a_s;
//   // arma::vec Omega = sample_Omega_cpp_noXb(Xbeta, n);
//   
//   arma::mat XtOmegaX = matrixProductXtOmegaX_SoR(X_arma, Y, X_centers, numTimeSpaceCov + ncov_psi, Omega,
//                                                  X_y_index, X_s_index);
//   
//   arma::vec knew = k - Omega % as_arma;
//   
//   arma::vec betaout = sample_beta_cpp_fast_sparse_SoR(X_arma, invB, b, knew, XtOmegaX,
//                                                       X_y_index, X_s_index, Y, X_centers, numTimeSpaceCov + ncov_psi);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));//,
//   // _["XtOmegaX"] = XtOmegaX));
// }
// 
// // [[Rcpp::export]]
// List sample_beta_omega_cpp(arma::vec beta,
//                            arma::mat &X, 
//                            arma::vec b, 
//                            arma::mat B, 
//                            arma::vec n, 
//                            arma::vec k){
//   
//   // sample Omega
//   arma::vec Omega = sample_Omega_cpp(X, beta, n);
//   
//   arma::mat tX = arma::trans(X);
//   arma::mat tXOmega = diagMatrixProd(tX, Omega);
//   arma::mat XtOmegaX = tXOmega * X;
//   
//   beta = sample_beta_cpp_fast(X, B, b, k, XtOmegaX);
//   
//   return(List::create(_["beta"] = beta,
//                       _["Omega"] = Omega));
// }
// 
// // [[Rcpp::export]]
// List sample_beta_omega_cpp_parallel(SEXP beta,
//                                     SEXP X, 
//                                     arma::vec b, 
//                                     arma::mat B, 
//                                     SEXP n, 
//                                     arma::vec k){
//   
//   NumericMatrix X_num = as<NumericMatrix>(X);
//   NumericVector beta_num = as<NumericVector>(beta);
//   NumericVector n_num = as<NumericVector>(n);
//   
//   // sample Omega
//   arma::vec Omega = sampleOmegaParallel(X_num, beta_num, n_num);
//   
//   arma::mat X_arma(X_num.begin(), X_num.nrow(), X_num.ncol(), false);
//   
//   arma::mat tX = arma::trans(X_arma);
//   arma::mat tXOmega = diagMatrixProd(tX, Omega);
//   arma::mat XtOmegaX = tXOmega * X_arma;
//   
//   arma::vec betaout = sample_beta_cpp_fast(X_arma, B, b, k, XtOmegaX);
//   
//   return(List::create(_["beta"] = betaout,
//                       _["Omega"] = Omega));
// }
// 
// // SAMPLER L (SCALE PARAMETER OF GAUSSIAN PROCESS)
// 
// // [[Rcpp::export]]
// double k_cpp(double x1, double x2, double a, double l){
//   // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
//   return a*exp(-(x1-x2)*(x1-x2)/(2*pow(l,2)));
//   // return 1;
// }
// 
// // [[Rcpp::export]]
// arma::mat K(arma::vec x1, arma::vec x2, double a, double l){
//   arma::mat res(x1.size(), x2.size());
//   
//   for(int i = 0; (unsigned)i < x1.size(); i++){
//     for(int j = 0; (unsigned)j < x2.size(); j++){
//       res(i,j) = k_cpp(x1[i],x2[j], a, l);
//     }  
//   }
//   
//   return res;
// }
// 
// // [[Rcpp::export]]
// double k2_cpp(arma::rowvec x1, arma::rowvec x2, double a, double l){
//   // return pow(1 + (x1-x2)*(x1-x2), - alphaGP);
//   return a*exp(-( pow(x1[0]-x2[0], 2) + pow(x1[1]-x2[1], 2) ) /(2*pow(l,2)));
// }
// 
// // [[Rcpp::export]]
// arma::mat K2(arma::mat x1, arma::mat x2, double a, double l){
//   arma::mat res(x1.n_rows, x2.n_rows);
//   
//   for(int i = 0; (unsigned)i < x1.n_rows; i++){
//     for(int j = 0; (unsigned)j < x2.n_rows; j++){
//       res(i,j) = k2_cpp(x1.row(i),x2.row(j), a, l);
//     }  
//   }
//   
//   return res;
// }
// 
// // [[Rcpp::export]]
// double loglikelihood_l_gp_cpp(double l, 
//                               double a, 
//                               int Y, 
//                               arma::mat XtOmegaX, 
//                               arma::mat Xtz){
//   
//   arma::vec years(Y);
//   std::iota(years.begin(), years.end(), 1);
//   arma::mat K_l = K(years, years, a, l);
//   
//   arma::mat invKl = arma::inv(K_l);
//   
//   arma::mat XtOmegaXpsolveK_l = XtOmegaX + invKl;
//   arma::mat identityMatrix = arma::eye(invKl.n_rows, invKl.n_rows);
//   arma::mat invKlXtOmegaX = invKl * XtOmegaX;
//   arma::mat IpXtOmegaXpsolveK_l = identityMatrix + invKlXtOmegaX;
//   
//   double logdet1det2;
//   double sign;
//   
//   log_det(logdet1det2, sign, IpXtOmegaXpsolveK_l);
//   
//   arma::vec term = arma::trans(Xtz) * arma::inv(XtOmegaXpsolveK_l) * Xtz;
//   
//   return (- .5 * logdet1det2 + .5 * term[0]);
// }
// 
// 
// // [[Rcpp::export]]
// arma::mat vec_subset_mat(const arma::mat& x, const arma::uvec& idx) {
//   return x.cols(idx);
// }
// 
// // [[Rcpp::export]]
// arma::vec subset_vec(const arma::vec& x, const arma::uvec& idx) {
//   return x.elem(idx);
// }
// 
// // [[Rcpp::export]]
// double sampler_l(double l, double a, arma::vec beta, arma::mat X,
//                  double a_l, double b_l, double sd_l,
//                  arma::vec a_s, int Y,
//                  arma::mat XtOmegaX, 
//                  arma::vec k, arma::vec Omega,
//                  arma::vec X_y_index){
//   
//   arma::vec idxes = arma::zeros(X.n_cols - Y);
//   idxes[0] = 0;
//   for(int i = 0; (unsigned)i < (idxes.size() - 1); i++) idxes[i + 1] = 1 + Y + i;
//   arma::uvec uidxes = arma::conv_to<arma::uvec>::from(idxes);
//   arma::mat X_no_y = vec_subset_mat(X, uidxes);
//   arma::vec beta_no_y = subset_vec(beta, uidxes);
//   
//   arma::vec c_i = X_no_y * beta_no_y + a_s;
//   
//   arma::vec y = k - Omega % c_i;
//   
//   arma::vec Xtz = arma::zeros(Y);
//   for(int i = 0; (unsigned)i < X_y_index.size(); i++){
//     Xtz[X_y_index[i] - 1] += y[i];
//   }
//   
//   double l_star = R::rnorm(l, sd_l);
//   
//   arma::mat XtOmegaX_subset = arma::zeros(Y, Y);
//   for(int i = 0; i < Y; i++){
//     XtOmegaX_subset(i,i) = XtOmegaX(i + 1, i + 1);
//   }
//   
//   if(l_star > 0){
//     
//     double loglikelihood_star = loglikelihood_l_gp_cpp(l_star, a, Y, XtOmegaX_subset, Xtz);
//     double loglikelihood_current = loglikelihood_l_gp_cpp(l, a, Y, XtOmegaX_subset, Xtz);
//     
//     double logprior_star = R::dgamma(l_star, a_l, 1 / b_l, 1);
//     double logprior_current = R::dgamma(l, a_l, 1 / b_l, 1);
//     
//     double logposterior_star = loglikelihood_star + logprior_star;
//     double logposterior_current = loglikelihood_current + logprior_current;
//     
//     if(R::runif(0, 1) < exp(logposterior_star - logposterior_current)){
//       
//       l = l_star;
//       
//     }
//     
//   }
//   
//   return l;
// }
// 
// // [[Rcpp::export]]
// double rcpp_log_dmvnorm_fast_cpp(arma::mat &inv_S, arma::vec &diag_S, 
//                                  double sigma_s, arma::vec &x) {
//   int n = x.size();
//   
//   double term1 = (-n / 2.0) * log(2 * M_PI) - sum(log(diag_S) + log(sigma_s));
//   
//   arma::vec term2 = inv_S * x;
//   
//   // arma::vec term3 = arma::trans(x) * term2;
//   
//   // double term4 = (.5) * (1 / (sigma_s*sigma_s)) * term3[0];
//   
//   return(0);
//   // return(term1 + term4);
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_l_grid_cpp(arma::vec l_s_grid, double sigma_s, 
//                             arma::cube &inv_K_s_grid, arma::mat &diag_K_s_grid,
//                             double a_l_S, double b_l_S, arma::vec a_s){
//   
//   arma::vec posterior_val = arma::zeros(l_s_grid.size());
//   
//   for(int j = 0; j < l_s_grid.size(); j++){
//     // for (j in 1:length(l_s_grid)) {
//     
//     double l_s = l_s_grid[j];
//     
//     arma::mat inv_K_s_grid_j = inv_K_s_grid.subcube(arma::span(), arma::span(), arma::span(j));
//     arma::vec diag_K_s_grid_j = arma::conv_to<arma::vec>::from(diag_K_s_grid.col(j));
//     // # K_s_grid_j <- K_s_grid[,,j] * sigma_s^2
//     // # inv_K_s_grid_j <- inv_K_s_grid[,,j] / sigma_s^2
//     // # diag_K_s_grid_j <- diag_K_s_grid[,j] * sigma_s
//     
//     double loglikelihood = 0;//rcpp_log_dmvnorm_fast_cpp(inv_K_s_grid_j, diag_K_s_grid_j, 
//     //                      sigma_s, a_s);
//     
//     // loglikelihood <- rcpp_log_dmvnorm_fast(inv_K_s_grid[,,j], 
//     // diag_K_s_grid[,j], sigma_s, a_s)
//     
//     // # loglikelihood <- rcpp_log_dmvnorm_fast(1, inv_K_s_grid_j, 
//     // #                                        diag_K_s_grid_j,  a_s)
//     //     
//     // # Sigma_l <- K2(X_tilde, X_tilde, sigma_s^2, l_s) + diag(exp(-10), nrow = nrow(X_tilde))
//     //     
//     // # (loglikelihood2 <- rcpp_log_dmvnorm( Sigma_l, rep(0, X_centers), a_s, F))
//     
//     double logPrior = R::dgamma(l_s, a_l_S, 1 / b_l_S, 1);
//     
//     posterior_val[j] = logPrior + loglikelihood;
//     
//   }
//   
//   return(posterior_val);
// }
// 
// // SAMPLE AS (STANDARD)
// 
// // [[Rcpp::export]]
// arma::vec sample_as_cpp_old(arma::vec k_s, arma::vec sites,
//                             arma::vec beta_psi, arma::mat &X_psi,
//                             arma::vec k,//arma::vec b_psi,//arma::vec c_psi, // 
//                             arma::vec z, arma::vec Omega,
//                             double sigma_a){
//   
//   arma::vec c_psi = X_psi * beta_psi;
//   arma::vec b_psi = k - Omega % c_psi;
//   
//   arma::vec a_s = arma::zeros(k_s.size());
//   
//   int index_site = 0;
//   
//   for(int i = 0; (unsigned)i < sites.size(); i++){
//     
//     int site = sites[i];
//     
//     int l = 1;
//     
//     // find rows associated with current site
//     if((unsigned)i != (sites.size() - 1)){
//       
//       while(k_s[index_site + l] == site){
//         l += 1;
//       }
//       
//     } else {
//       
//       l = k_s.size() - index_site;
//       
//     }
//     
//     IntegerVector indexes_site(l);
//     for(int j = 0; j < l; j++){
//       indexes_site[j] = index_site + j;
//     }
//     index_site += l;
//     
//     // arma::mat data_psi_site(l, X_psi.n_cols);
//     // arma::vec z_site = arma::zeros(l);
//     arma::vec w_site = arma::zeros(l);
//     // arma::vec c_site = arma::zeros(l);
//     arma::vec b_site = arma::zeros(l);
//     for(int j = 0; j < l; j++){
//       // data_psi_site.row(j) = X_psi.row(indexes_site[j]);
//       // z_site[j] = z[indexes_site[j]];
//       w_site[j] = Omega[indexes_site[j]];
//       // c_site[j] = c_psi[indexes_site[j]];
//       b_site[j] = b_psi[indexes_site[j]];
//     }
//     
//     // sample x
//     // arma::vec k_site = (z_site - .5);
//     // arma::vec n_site = arma::ones(l);
//     // arma::vec c_site = data_psi_site * beta_psi;
//     
//     double a = 1 / (sum(w_site) + 1 / (sigma_a*sigma_a));
//     // double b = sum(k_site - w_site % c_site);
//     double b = sum(b_site);
//     
//     double x = R::rnorm(b * a, sqrt(a));
//     
//     for(int j = 0; j < l; j++){
//       a_s[indexes_site[j]] = x;
//     }
//     
//   }
//   
//   return(a_s);
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_eps_cpp(arma::vec k_s, arma::vec sites,
//                          arma::vec &c_psi,
//                          arma::vec k,//arma::vec b_psi,//arma::vec c_psi, // 
//                          arma::vec z, arma::vec Omega,
//                          double sigma_a){
//   
//   arma::vec b_psi = k - Omega % c_psi;
//   
//   arma::vec a_s = arma::zeros(k_s.size());
//   
//   int index_site = 0;
//   
//   for(int i = 0; (unsigned)i < sites.size(); i++){
//     
//     int site = sites[i];
//     
//     int l = 1;
//     
//     // find rows associated with current site
//     if((unsigned)i != (sites.size() - 1)){
//       
//       while(k_s[index_site + l] == site){
//         l += 1;
//       }
//       
//     } else {
//       
//       l = k_s.size() - index_site;
//       
//     }
//     
//     IntegerVector indexes_site(l);
//     for(int j = 0; j < l; j++){
//       indexes_site[j] = index_site + j;
//     }
//     index_site += l;
//     
//     // arma::mat data_psi_site(l, X_psi.n_cols);
//     // arma::vec z_site = arma::zeros(l);
//     arma::vec w_site = arma::zeros(l);
//     // arma::vec c_site = arma::zeros(l);
//     arma::vec b_site = arma::zeros(l);
//     for(int j = 0; j < l; j++){
//       // data_psi_site.row(j) = X_psi.row(indexes_site[j]);
//       // z_site[j] = z[indexes_site[j]];
//       w_site[j] = Omega[indexes_site[j]];
//       // c_site[j] = c_psi[indexes_site[j]];
//       b_site[j] = b_psi[indexes_site[j]];
//     }
//     
//     // sample x
//     // arma::vec k_site = (z_site - .5);
//     // arma::vec n_site = arma::ones(l);
//     // arma::vec c_site = data_psi_site * beta_psi;
//     
//     double a = 1 / (sum(w_site) + 1 / (sigma_a*sigma_a));
//     // double b = sum(k_site - w_site % c_site);
//     double b = sum(b_site);
//     
//     double x = R::rnorm(b * a, sqrt(a));
//     
//     for(int j = 0; j < l; j++){
//       a_s[indexes_site[j]] = x;
//     }
//     
//   }
//   
//   return(a_s);
// }
// 
// double prior_as(double a_s, double sum_as, double beta, double a){
//   return(pow(a_s*a_s + sum_as + beta,-a));
// }
// 
// // [[Rcpp::export]]
// arma::vec sample_as_mh_t_cpp(arma::vec a_s, arma::vec k_s, arma::vec sites, arma::vec Xbeta,
//                              arma::vec z, double a_sigma, double b_sigma){
//   
//   
//   int index_site = 0;
//   
//   double sum_as = 0;
//   
//   for(int i = 0; (unsigned)i < sites.size(); i++){
//     
//     int l = 1;
//     
//     int site = sites[i];
//     
//     // find rows associated with current site
//     if((unsigned)i != (sites.size() - 1)){
//       
//       while(k_s[index_site + l] == site){
//         l += 1;
//       }
//       
//     } else {
//       
//       l = k_s.size() - index_site;
//       
//     }
//     
//     IntegerVector indexes_site(l);
//     for(int j = 0; j < l; j++){
//       indexes_site[j] = index_site + j;
//     }
//     index_site += l;  
//     
//     sum_as += a_s[indexes_site[0]] * a_s[indexes_site[0]];
//     
//   }
//   
//   index_site = 0;
//   
//   for(int i = 0; (unsigned)i < sites.size(); i++){
//     
//     int site = sites[i];
//     
//     int l = 1;
//     
//     // find rows associated with current site
//     if((unsigned)i != (sites.size() - 1)){
//       
//       while(k_s[index_site + l] == site){
//         l += 1;
//       }
//       
//     } else {
//       
//       l = k_s.size() - index_site;
//       
//     }
//     
//     IntegerVector indexes_site(l);
//     for(int j = 0; j < l; j++){
//       indexes_site[j] = index_site + j;
//     }
//     index_site += l;  
//     
//     arma::vec c_site = arma::zeros(l); 
//     arma::vec z_site = arma::zeros(l);
//     for(int j = 0; j < l; j++){
//       c_site[j] = Xbeta[indexes_site[j]];
//       z_site[j] = z[indexes_site[j]];
//     }  
//     
//     // sample x
//     arma::vec k_site = (z_site - .5);
//     arma::vec n_site = arma::ones(l);
//     
//     double a_s_old = a_s[indexes_site[0]];
//     double a_s_new = R::rnorm(a_s_old, .05);
//     
//     double loglikelihood_current = 0;
//     for(int j = 0; j < l; j++){
//       loglikelihood_current += z[indexes_site[j]] * (Xbeta[indexes_site[j]] + a_s_old) - 
//         log(1 + exp(Xbeta[indexes_site[j]] + a_s_old));
//     }
//     
//     double loglikelihood_new = 0;
//     for(int j = 0; j < l; j++){
//       loglikelihood_new += z[indexes_site[j]] * (Xbeta[indexes_site[j]] + a_s_new) - 
//         log(1 + exp(Xbeta[indexes_site[j]] + a_s_new));
//     }
//     
//     sum_as -= (a_s_old * a_s_old);
//     
//     double prior_current = prior_as(a_s_old, sum_as, b_sigma, (a_sigma + sites.size())/ 2);
//     double prior_new = prior_as(a_s_new, sum_as, b_sigma, (a_sigma + sites.size())/ 2);
//     
//     double logposterior_current = dt_cpp(a_s_old, 2 * a_sigma, 0, b_sigma / a_sigma, 1) +
//       loglikelihood_current;
//     double logposterior_new = dt_cpp(a_s_new, 2 * a_sigma, 0, b_sigma / a_sigma, 1) +
//       loglikelihood_new;
//     
//     if(R::runif(0,1) < exp(logposterior_new - logposterior_current)){
//       a_s_old = a_s_new;
//     }
//     
//     sum_as += (a_s_old * a_s_old);
//     
//     for(int j = 0; j < l; j++){
//       a_s[indexes_site[j]] = a_s_old;
//     }
//     
//     // mu_clusters[site - 1] = b * a;
//     // sigma_clusters[site - 1] = sqrt(a);
//     
//   }
//   
//   return(a_s);
// }
// 
// // CREATE GRID
// 
// // [[Rcpp::export]]
// bool isPointInBandRight(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
//   
//   for(int k = 0; k < X_tilde.n_rows; k++){
//     
//     if((X_tilde(k,1) < y_grid[j + 1]) & (X_tilde(k,1) > y_grid[j - 1])){
//       if(X_tilde(k,0) < x_grid[i + 1]){
//         return(true);
//       }
//     } 
//     
//   }
//   
//   return(false);
// }
// 
// // [[Rcpp::export]]
// bool isPointInBandLeft(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j) {
//   
//   for(int k = 0; k < X_tilde.n_rows; k++){
//     
//     if((X_tilde(k,1) < y_grid[j + 1]) & (X_tilde(k,1) > y_grid[j - 1])){
//       if(X_tilde(k,0) > x_grid[i - 1]){
//         return(true);
//       }
//     } 
//     
//   }
//   
//   return(false);
// }
// 
// // [[Rcpp::export]]
// bool isPointInBandUp(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
//   
//   for(int k = 0; k < X_tilde.n_rows; k++){
//     
//     if((X_tilde(k,0) < x_grid[i + 1]) & (X_tilde(k,0) > x_grid[i - 1])){
//       if(X_tilde(k,1) > y_grid[j-1]){
//         return(true);
//       }
//     }
//     
//   }
//   
//   return(false);
//   
// }
// 
// // [[Rcpp::export]]
// bool isPointInBandDown(arma::mat X_tilde, arma::vec x_grid, arma::vec y_grid, int i, int j){
//   
//   for(int k = 0; k < X_tilde.n_rows; k++){
//     
//     if((X_tilde(k,0) < x_grid[i + 1]) & (X_tilde(k,0) > x_grid[i - 1])){
//       if(X_tilde(k,1) < y_grid[j+1]){
//         return(true);
//       }
//     }
//     
//   }
//   
//   return(false);
//   
// }
// 
// // [[Rcpp::export]]
// IntegerVector findClosestPoint(arma::mat XY_sp, arma::mat X_tilde){
//   
//   IntegerVector closestPoint(XY_sp.n_rows);
//   
//   for(int k = 0; k < XY_sp.n_rows; k++){
//     
//     double newDistance = 0;
//     double minDistance = exp(10);
//     int bestIndex = 0;
//     
//     for(int i = 0; i < X_tilde.n_rows; i++){
//       newDistance = pow(X_tilde(i, 0) - XY_sp(k, 0), 2) + pow(X_tilde(i, 1) - XY_sp(k, 1), 2);
//       
//       if(newDistance < minDistance){
//         minDistance = newDistance;
//         bestIndex = i + 1;
//       }
//     }
//     
//     closestPoint[k] = bestIndex;
//     
//   }
//   
//   return(closestPoint);
// }
// 
// //// HYPERPARAMETER SAMPLING FUNCTIONS
// 
// 
// // [[Rcpp::export]]
// arma::mat matrixProductXtOmegaX_year(int Y, arma::vec Omega,
//                                      arma::vec X_y_index){
//   
//   arma::mat XtOmegaX2 = arma::zeros(Y, Y);
//   
//   // year covariates times year covariates
//   for (int i = 0; (unsigned)i < X_y_index.size(); i++){
//     
//     XtOmegaX2(X_y_index[i] - 1, X_y_index[i] - 1) += Omega[i];
//     
//   }
//   
//   return(XtOmegaX2);
// }
// 
// // [[Rcpp::export]]
// arma::mat matrixProductXtOmegaX_spatial(int X_centers, arma::vec Omega,
//                                         arma::vec X_s_index){
//   
//   arma::mat XtOmegaX2 = arma::zeros(X_centers, X_centers);
//   
//   // spatial  covariates times spatial covariates
//   for (int i = 0; (unsigned)i < X_s_index.size(); i++){
//     
//     XtOmegaX2(X_s_index[i] - 1, X_s_index[i] - 1) += Omega[i];
//     
//   }
//   
//   return(XtOmegaX2);
// }
// 
// // [[Rcpp::export]]
// arma::vec XpsiYz(arma::vec X_y_index, arma::vec z, int Y){
//   
//   arma::vec out = arma::zeros(Y);
//   
//   for(int i = 0; i < X_y_index.size(); i++){
//     out(X_y_index[i] - 1) += z[i];
//   }
//   
//   return(out);
// }
// 
// // [[Rcpp::export]]
// arma::vec XpsinoYbetaz(arma::vec X_s_index, int X_centers, int ncov,
//                        arma::mat &X_cov, arma::vec& beta){
//   
//   // X_psi[,-(1:Y),drop=F] %*% beta_psi[-(1:Y)]
//   
//   return(X_cov * beta);
//   // arma::vec out = arma::zeros(X.n_cols);
//   
//   // for(int i = 0; i < X_s_index.size(); i++){
//   //   Xk(Y + X_s_index[i] - 1) += k[i];
//   // }
//   // for(int i = 0; i < ncov; i++){
//   //   Xk(Y + centers + i) = as_scalar(k.t() * X.col(Y + centers + i));
//   // }
//   
// }
// 
// 
// //// COMPUTE LIKELIHOOD
// 
// // [[Rcpp::export]]
// arma::vec computelikelihood_cpp(arma::vec Occs, arma::vec p, arma::vec z_all){
//   
//   arma::vec likelihoods = arma::zeros(Occs.size());
//   
//   for(int i = 0; i < Occs.size(); i++){
//     if(z_all[i] == 0){
//       likelihoods[i] = 1;
//     } else {
//       likelihoods[i] = R::dbinom(Occs[i], 1, p[i], 0);
//     }
//   }
//   
//   return(likelihoods);
// }

