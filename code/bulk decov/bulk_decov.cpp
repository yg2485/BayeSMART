#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <RcppDist.h>

// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;

// static int num(IntegerVector x, int c);
Rcpp::List bulk_decov(arma::mat X, arma::mat Y, IntegerVector cell_type, 
    IntegerVector gamma, NumericVector v,  NumericVector s, NumericVector g,  
    int iter, int burn, double tau_pi, double tau_mu, int d);

// set seed
// [[Rcpp::export]]
void set_seed(double seed) {
  Rcpp::Environment base_env("package:base");
  Rcpp::Function set_seed_r = base_env["set.seed"];
  set_seed_r(std::floor(std::fabs(seed)));
}


// [[Rcpp::export]]
Rcpp::List imoscato_new(
    arma::mat X, 
    arma::mat Y, 
    IntegerVector cell_type, 
    IntegerVector gamma, 
    NumericVector v,  
    NumericVector s, 
    NumericVector g, 
    int iter, 
    int burn,
    double tau_pi = 0.01,
    double tau_mu = 0.1,
    int d = 42
){
  set_seed(d);
  
  // Read data information
  int n = Y.n_rows;
  int p = Y.n_cols;
  int m = X.n_rows;
  int K = Rcpp::unique(cell_type).size(); // here K is the spatial domain size
  Rcout<<K<< " is the number of domains\n";
  
  NumericVector number(K);
  NumericVector prop(K);
  for (int k = 0; k < K; k++){
    // number(k) = num(cell_type, k+1);
    // prop(k) = number(k)/m;
    prop(k) = pow(K, -1);
  }

  
  // Hyperparamters
  double a_mu = 0.001;
  double b_mu = 0.001;
  double a_0 = 0.001;
  double b_0 = 0.001;
  
  
  // Set temporary variables

  int t, i, j, k, kk, num = 10;
  double accept_mu = 0, h_temp, accept_pi = 0;
  double try_mu = 0, try_pi = 0;
  double hastings, mu_sum_temp, mu_sum_new;


  NumericMatrix M_new(K, p);
  
  NumericVector pi_temp(K);
  NumericVector pi_new(K);

  
  // Create the spaces to store the results
  arma::cube Pi_store(n, K, iter);
  arma::cube M_store(K, p, iter);

  
  // Initialization
  IntegerVector A(K);
  for (k = 0; k < K; k++){
    A(k) = 1;
  }
  
  NumericMatrix ct_p(n, K);
  for (k = 0; k < K; k++){
    for (i = 0; i < n; i++){
      ct_p(i, k) = prop(k);
      // ct_p(i, k) = pow(K, -1);
    }
  }
  
  
  NumericMatrix M(K, p);
  for (k = 0; k < K; k++){
    for (j = 0; j < p; j++){
      M(k, j) = mean(mean(X)); //Initialize of mu_kj
    }
  }
  Rcout<< mean(mean(X)) << " is the value of the initial mu\n";
  
  
  IntegerVector clusters(K-1);
  for (k = 0; k < (K-1); k++){
    clusters(k) = k;
  } // this is for updating Pi
  
  
  // Start MCMC algorithms
  for (t = 0; t < iter; t++){

    
    // RWMH to update the cell type proportion Pi ##############################
    for (i = 0; i < n; i++) {
      // Propose new pi vector
      for (kk = 0; kk < K; kk++){
        pi_temp(kk) = ct_p(i, kk);
        pi_new(kk) = ct_p(i, kk);
      }

      int k_temp = Rcpp::as<int>(Rcpp::RcppArmadillo::sample(clusters, 1, TRUE));
      h_temp = r_truncnorm(0, tau_pi, -pi_temp(k_temp + 1),  pi_temp(k_temp));
      // Rcout<<h_temp<< "h\n";

      pi_new(k_temp) = pi_temp(k_temp) - h_temp;
      pi_new(k_temp + 1) = pi_temp(k_temp + 1) + h_temp;


      // Calculate hastings ratio
      hastings = 0;

      // SRT data likelihood ratio
      for (j = 0; j < p; j++){
        if (gamma(j) == 1){
          mu_sum_new = sum(pi_new*M(_, j));
          mu_sum_temp = sum(pi_temp*M(_, j));
          // Rcout<<mu_sum_new<< "mu_sum_new\n";
          // Rcout<<mu_sum_temp<< "mu_sum_temp\n";

          hastings = hastings +  Y(i, j)*(log(mu_sum_new*s(i)*g(j))) - mu_sum_new*s(i)*g(j);
          hastings = hastings - Y(i, j)*(log(mu_sum_temp*s(i)*g(j))) + mu_sum_temp*s(i)*g(j);

        }
      }

      // Prior ratio
      // for (k = 0; k < K; k++) {
      //
      //   hastings = hastings + (A(k) - 1)*log(pi_new(k));
      //   hastings = hastings - (A(k) - 1)*log(pi_temp(k));
      // }

      // Rcout<<hastings<< "is the ratio\n";

      // check if accept the proposed new values
      if (hastings > log(unif_rand())){
        for (k = 0; k < K; k++){
          ct_p(i, k) = pi_new(k);
        }
        if (t >= burn){
          accept_pi = accept_pi + 1;
          }
      }
      if (t >= burn){
        try_pi = try_pi + 1;
        }

    } // End of updating the cell type proportion Pi #########################





    // RWMH to update the normalized expression levels M #######################
    for (j = 0; j < p; j++){
      if (gamma(j) == 0) {
        double mu_temp = M(0, j);
        double mu_new = exp(rnorm(1, log(mu_temp),  tau_mu)(0));

        hastings = 0;

        // SRT data likelihood ratio
        for (i = 0; i < n; i++) {
          mu_sum_temp = mu_temp;
          mu_sum_new =  mu_new;

          hastings = hastings +  Y(i, j)*(log(mu_sum_new*s(i)*g(j))) - mu_sum_new*s(i)*g(j);
          hastings = hastings - Y(i, j)*(log(mu_sum_temp*s(i)*g(j))) + mu_sum_temp*s(i)*g(j);
        }
        

        //scRNA-seq data likelihood ratio
        for (i = 0; i < m; i++) {
          hastings = hastings + X(i, j)*(log(mu_new*v(i))) - mu_new*v(i);
          hastings = hastings - X(i, j)*(log(mu_temp*v(i))) + mu_temp*v(i);
        }

        // Prior ratio
        hastings = hastings + (a_0 - 1)*log(mu_new) - b_0*mu_new;
        hastings = hastings - ((a_0 - 1)*log(mu_temp) - b_0*mu_temp);

        // check if accept the proposed new values
        if (hastings > log(unif_rand())) {
          for (k = 0; k < K; k++) {
            M(k, j) = mu_new;
          }

          if (t >= burn){
            accept_mu = accept_mu + 1;
          }
        }

        if (t >= burn){
          try_mu = try_mu + 1;
        }
      } // End of updating M for non-discriminating genes

      else{
        for (k = 0; k < K; k++) {
          double mu_temp = M(k, j);
          double mu_new = exp(rnorm(1, log(mu_temp), tau_mu)(0));

          M_new(_, j) = M(_, j);
          M_new(k, j) = mu_new;

          hastings = 0;

          // SRT data likelihood ratio
          for (i = 0; i < n; i++) {
            mu_sum_temp =  sum(ct_p(i, _)*M(_, j));
            mu_sum_new =  sum(ct_p(i, _)*M_new(_, j));

            hastings = hastings +  Y(i, j)*(log(mu_sum_new*s(i)*g(j))) - mu_sum_new*s(i)*g(j);
            hastings = hastings - Y(i, j)*(log(mu_sum_temp*s(i)*g(j))) + mu_sum_temp*s(i)*g(j);
          }


          // scRNA-seq data likelihood ratio
          for (i = 0; i < m; i++) {
            if (cell_type(i) == k + 1) {
              hastings = hastings + X(i, j)*(log(mu_new*v(i))) - mu_new*v(i);
              hastings = hastings - X(i, j)*(log(mu_temp*v(i))) + mu_temp*v(i);
            }
          }

          // Prior ratio
          hastings = hastings + (a_mu - 1)*log(mu_new) - b_mu*mu_new;
          hastings = hastings - ((a_mu - 1)*log(mu_temp) - b_mu*mu_temp);

          // Check if accept the proposed new values
          if (hastings > log(unif_rand())) {
            M(k, j) = mu_new;

            if (t >= burn){
              accept_mu = accept_mu + 1;
            }
          }

          if (t >= burn){
            try_mu = try_mu + 1;
          }
        }

      } // End of updating M for discriminating genes
    } // End of updating the normalized expression levels M ######################



    // Monitor the process
    if((t*100/(iter-1)) == num) {
      Rcout<<num<< "% has been done\n";
      num = num + 10;
    }
    
    
    // Store the results
    Pi_store.slice(t) = as<arma::mat>(ct_p);
    M_store.slice(t) = as<arma::mat>(M);
    
    
  }// End of iterations
  

  accept_mu = accept_mu/try_mu;
  accept_pi = accept_pi/try_pi;

  return Rcpp::List::create(Rcpp::Named("Pi_store") = Pi_store,
                            Rcpp::Named("M_store") = M_store,
                            Rcpp::Named("accept_mu") = accept_mu,
                            Rcpp::Named("accept_pi") = accept_pi);
  
} // End of function


