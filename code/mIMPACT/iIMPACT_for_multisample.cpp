#include <RcppArmadillo.h>
#include <RcppDist.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo,RcppDist)]]
using namespace Rcpp;
using namespace arma;


/*
 define the type of function in advance
 */
arma::mat rdirichlet_cpp(NumericVector alpha, int size);
double lgamma_function(double alpha);
Rcpp::List mIMPACT(arma::mat V, arma::mat Y, arma::mat P,  int K, NumericVector e, double f, 
                   NumericVector eta, double tau, double alpha_gamma, NumericMatrix beta_gamma, NumericVector mu_initial, 
                   NumericVector sigma_initial, NumericVector alpha, NumericVector omega_initial, arma::mat Z_initial, NumericVector weight);
arma::mat rnorm_function(int n, NumericVector mu, arma::mat sigma);
double dnorm_function(arma::mat alpha, NumericVector mu, arma::mat sigma);
double lbeta_function(arma::mat alpha);
// arma::vec rMultivariateNormal(arma::vec mean, arma::mat cov);
// arma::mat rInvWishart(double df, arma::mat Sigma);
arma::mat rinvwish(double v, arma::mat S);
arma::vec mvrnorm(vec mu, mat Sigma);


// [[Rcpp::export]]
Rcpp::List mIMPACT(arma::mat V, arma::mat Y, arma::mat P,  int K, NumericVector e, NumericVector f, 
                   NumericVector eta, double tau, double alpha_gamma, NumericMatrix beta_gamma, 
                   NumericVector mu_initial, NumericMatrix sigma_initial, NumericVector alpha, 
                   NumericVector omega_initial, arma::mat Z_initial, NumericVector weight) {
  // Read data information
  int N = P.n_rows;
  int n_neighbor = P.n_cols;
  int Y_dim = Y.n_cols;
  int V_dim = V.n_cols;
  
  
  // Set algorithm settings
  int iter = 20000;
  //int burn = iter/2;
  
  int q, qq, qqq, h, count;
  int count_2 = 10;
  
  NumericMatrix mu_store(iter, Y_dim*K);
  NumericMatrix mu(K, Y_dim);
  arma::cube sigma(K, Y_dim, Y_dim); // revise
  NumericMatrix sigma_store(iter, K*Y_dim*Y_dim); // revise
  NumericMatrix Z_store(iter, N);
  NumericVector Z(N);
  NumericMatrix omega_store(iter, V_dim*K);
  NumericMatrix omega(K, V_dim);
  
  NumericVector count_temp(K);
  arma::mat average_temp(K, Y_dim);
  
  NumericVector V_temp(V_dim);
  arma::mat omega_temp(1, V_dim);
  NumericVector sample_index(K);
  
  // Initialization
  for(q = 0; q < K; q++)
  {
    for (qq = 0; qq < Y_dim; qq++)
    {
      mu(q, qq) = mu_initial(qq);
      for(qqq = 0; qqq < Y_dim; qqq++){
        sigma(q, qq, qqq) = sigma_initial(qq, qqq);
      }
      
    }
    sample_index(q) = q + 1; 
  }
  
  for (q = 0; q < N; q++)
  {
    Z(q) = Z_initial(q, 0);
  }
  
  for(q = 0; q < K; q++)
  {
    for (qq = 0; qq < V_dim; qq++)
    {
      omega(q, qq) = omega_initial(qq);
    }
    sample_index(q) = q + 1; 
  }
  
  
  
  // MCMC
  for(int t = 0; t < iter; t++)
  {
    if (t*100/iter == count_2)
    {
      Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }

    // reset to zero
    for(q = 0; q < K; q++)
    {
      count_temp(q) = 0;
      for (qq = 0; qq < Y_dim; qq++){
        average_temp(q, qq) = 0;
      }
    }

    // compute count and y_average
    for(q = 0; q < K; q++)
    {
      for (h = 0; h < N; h++){
        if (Z(h) == q + 1){
          count_temp(q) ++;
          for (qq = 0; qq < Y_dim; qq++){
            average_temp(q, qq) += Y(h, qq); // get sum y_k
          }
        }
      }
    }

    for(q = 0; q < K; q++)
    {
      for (qq = 0; qq < Y_dim; qq++){
        if (count_temp(q) > 0){
          average_temp(q, qq) =  average_temp(q, qq)/count_temp(q);} //get y_k bar
      }
    }


    // Update sigma
    for(q = 0; q < K; q++){
      double alpha_gamma_new = alpha_gamma + count_temp(q);
      arma::mat beta_gamma_new(Y_dim, Y_dim);

      for (qq = 0; qq < Y_dim; qq++){
        for (qqq = 0; qqq < Y_dim; qqq++){
          double sum_temp3 = beta_gamma(qq,qqq);
          for (h = 0; h < N; h++){
            if (Z(h) == q + 1){
              sum_temp3 += (Y(h, qq) - average_temp(q, qq))*(Y(h, qqq) - average_temp(q, qqq));
            }
          }
          sum_temp3 += tau*count_temp(q)*(average_temp(q, qq) - eta(qq))*(average_temp(q, qqq) - eta(qqq))/(tau + count_temp(q));
          beta_gamma_new(qq, qqq) = sum_temp3;
        }
      }
      arma::mat new_sigma = rinvwish(alpha_gamma_new, beta_gamma_new);
      for (qq = 0; qq < Y_dim; qq++){
        for (qqq = 0; qqq < Y_dim; qqq++){
          sigma(q, qq, qqq) = new_sigma(qq, qqq);
        }
      }
    }


    // Update mu
    for(q = 0; q < K; q++)
    {
      arma::vec eta_new(Y_dim);
      for(qq = 0; qq < Y_dim; qq++)
      {
        eta_new(qq) = (eta(qq)*tau + count_temp(q)*average_temp(q, qq))/(count_temp(q) + tau);
      }
      arma::mat sigma_new(Y_dim, Y_dim);
      for(int i=0; i < Y_dim; i++){
        for(int j=0; j < Y_dim; j++){
          sigma_new(i,j) = sigma(q, i, j)/(tau + count_temp(q));
        }
      }
      arma::vec mu_new = mvrnorm(eta_new, sigma_new);
      for(qq = 0; qq < Y_dim; qq++){
        mu(q,qq) = mu_new(qq);
      }
    }


    // Update omega
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < V_dim; qq++)
      {
        V_temp(qq) = alpha(qq);
      }

      for(qq = 0; qq < N; qq++)
      {
        if (Z(qq) == q + 1){
          for(qqq = 0; qqq < V_dim; qqq++){
            V_temp(qqq) = V_temp(qqq) + V(qq, qqq);
          }}
      }
      omega_temp = rdirichlet_cpp(V_temp, V_dim);

      for(qq = 0; qq < V_dim; qq++)
      {
        omega(q, qq) = omega_temp(0, qq);
      }
    }


    // Update Z
    for (h = 0; h < N; h++){
      double sum_temp = 0;
      double sum_temp2 = 0;

      NumericVector prob_temp(K);
      NumericVector neighbor_index(K);

      for (q = 0; q < n_neighbor; q++){
        if (P(h, q) != 0){
          neighbor_index(Z(P(h, q) - 1) - 1) = neighbor_index(Z(P(h, q) - 1) - 1) + 1;
        }
      }

      for(q = 0; q < K; q++)
      {
        arma::mat new_sigma(Y_dim, Y_dim);
        for (qq = 0; qq < Y_dim; qq++){
          for (qqq = 0; qqq < Y_dim; qqq++){
            new_sigma(qq, qqq) = sigma(q, qq, qqq);
          }
        }
        arma::mat sigma_inv = arma::inv(new_sigma); // get the inverse of sigma

        prob_temp(q) = e(q) + f(h)*neighbor_index(q);

        for (qq = 0; qq < Y_dim; qq++){
          double record = 0;
          for (qqq= 0; qqq < Y_dim; qqq++){
            record += (Y(h, qqq) - mu(q, qqq))*sigma_inv(qq,qqq);
          }
          prob_temp(q) = prob_temp(q) - 0.5*record*(Y(h, qq) - mu(q, qq));
        }

        prob_temp(q) -= 0.5*log(arma::det(new_sigma));

        for (qq = 0; qq < V_dim; qq++){
          prob_temp(q) = prob_temp(q) + weight(h)*V(h, qq)*log(omega(q, qq) + 0.00000000000000000001);
        }
        sum_temp = sum_temp + prob_temp(q);
      }

      for(q = 0; q < K; q++)
      {
        prob_temp(q) = prob_temp(q) - sum_temp/K;
        if (prob_temp(q) > 709){prob_temp(q) = 709;}
        if (prob_temp(q) < -300){prob_temp(q) = -300; }
        prob_temp(q) = exp(prob_temp(q));
        sum_temp2 = sum_temp2 + prob_temp(q);
      }

      if(sum_temp2 == 0){
        for(q = 0; q < K; q++)
        {
          prob_temp(q) = 1.0/K;
        }
        Rcout <<t<< "% zero\n";
      }
      else{
        for(q = 0; q < K; q++)
        {
          prob_temp(q) = prob_temp(q)/sum_temp2 + 0.0000000000000000001;
        }
      }

      Z(h) = sample(sample_index, 1, 1==0, prob_temp)(0);
    }


    // Monitor the process
    if (t*100/iter == count_2)
    {
      Rcout <<count_2<< "% has been done\n";
      count_2 = count_2 + 10;
    }
    count = 0;
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < Y_dim; qq++)
      {
        mu_store(t, count) = mu(q, qq);
        count++;
      }
    }

    count= 0;
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < Y_dim; qq++)
      {
        for (qqq = 0; qqq < Y_dim; qqq++){
          sigma_store(t, count) = sigma(q, qq, qqq); // record by rows
          count++;
        }
      }
    }

    count = 0;
    for(q = 0; q < K; q++)
    {
      for(qq = 0; qq < V_dim; qq++)
      {
        omega_store(t, count) = omega(q, qq);
        count++;
      }
    }

    for(q = 0; q < N; q++)
    {
      Z_store(t, q) = Z(q);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("mu") = mu_store, Rcpp::Named("sigma") = sigma_store, Rcpp::Named("omega") = omega_store, Rcpp::Named("Z") = Z_store);
}



// [[Rcpp::export]]
arma::mat rdirichlet_cpp(NumericVector alpha, int size) {
  
  arma::mat sample = arma::zeros(1, size);
  
  double sum_temp = 0;
  for (int j = 0; j < size; ++j){
    double cur = rgamma(1, alpha(j), 1.0)(0);
    sample(0, j) = cur; 
    sum_temp += cur; 
  }
  
  for (int j = 0; j < size; ++j){
    sample(0, j) = sample(0, j)/sum_temp; 
  }
  return(sample);
}

// [[Rcpp::export]]
double lbeta_function(arma::mat alpha) {
  double re = 0.0;
  double sum_temp = 0.0;
  
  int N =alpha.n_cols;
  for (int i = 0; i < N; i++){
    re += lgamma_function(exp(alpha(0, i)));
    sum_temp += exp(alpha(0, i));
  }
  
  re -= lgamma_function(sum_temp);
  
  return(re);
}


// [[Rcpp::export]]
double lgamma_function(double alpha) {
  NumericVector alpha_vec(1);
  
  alpha_vec(0) = alpha;
  double qq = Rcpp::sum(Rcpp::lgamma(alpha_vec));
  
  return(qq);
}


// [[Rcpp::export]]
double dnorm_function(arma::mat alpha, NumericVector mu, arma::mat sigma) {
  arma::vec qq = dmvnorm(alpha, mu, sigma);
  
  return(qq(0, 0));
}


// [[Rcpp::export]]
arma::mat rnorm_function(int n, NumericVector mu, arma::mat sigma) {
  arma::mat qq = rmvnorm(n, mu, sigma);
  return(qq);
}

// the following codes are from 
// https://github.com/fditraglia/econ722/blob/master/RcppArmadillo/InverseWishartSampler/InverseWishart.cpp


// [[Rcpp::export]]
arma::mat rinvwish(double v, arma::mat S){
  RNGScope scope;
  int p = S.n_rows;
  arma::mat L = chol(inv_sympd(S), "lower");
  arma::mat sims(p, p, arma::fill::zeros);
  arma::mat A(p, p, arma::fill::zeros);
  for(int i = 0; i < p; i++){
    int df = v - (i + 1) + 1; //zero-indexing
    A(i,i) = sqrt(R::rchisq(df)); 
  }
  for(int row = 1; row < p; row++){
    for(int col = 0; col < row; col++){
      A(row, col) = R::rnorm(0,1);
    }
  }
  arma::mat LA_inv = inv(trimatl(trimatl(L) * trimatl(A)));
  sims = LA_inv.t() * LA_inv;
  return(sims);
}


// [[Rcpp::export]]
vec mvrnorm(vec mu, mat Sigma) {
  RNGScope scope;
  int p = Sigma.n_cols;
  vec x = rnorm(p);
  vec eigval;
  mat eigvec;
  eig_sym(eigval, eigvec, Sigma);
  x = eigvec * diagmat(sqrt(eigval)) * x;
  x += mu;
  return x;
}

