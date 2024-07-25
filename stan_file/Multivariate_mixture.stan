data {
  int N; //number of observations
  int K; //number of classes
  int NY; //number of observed variables
  array[N] vector[NY] Y; //data
}

parameters {
  simplex[K] theta; //mixing proportions
  array[NY] positive_ordered[K] mu;; //mixture component means
  //cholesky_factor_corr[NY] L[K]; //cholesky factor of covariance
  array[K] corr_matrix[NY] Rho; // Correlation matix
  array[K] vector<lower=0>[NY] sigma;
  
  positive_ordered[K] p;
}
transformed parameters {
  array[K] vector[NY] mu_ord;
  for (j in 1:NY) mu_ord[:, j] = to_array_1d(mu[j, :]);
}
model {
  p ~ normal(0,1);
  theta ~ beta(1, 1);
  
  for(k in 1:K){
    mu[, k] ~ normal(p[k], 2);
    // L[k] ~ lkj_corr_cholesky(5);
    Rho[k] ~ lkj_corr(5);
    sigma[k] ~ std_normal();
  }

 {
    vector[K] ps;
    for (n in 1:N){
      for (k in 1:K){
        ps[k] = log(theta[k]) + multi_normal_lpdf(Y[n] | to_vector(mu_ord[k]), quad_form_diag(Rho[k], sigma[k]));
//         ps[k] = log(theta[k]) + multi_normal_cholesky_lpdf(Y[n] | // to_vector(mu_ord[k]), diag_pre_multiply(sigma[k], L[k]));

//increment log probability of the gaussian
      }
      target += log_sum_exp(ps);
    }
  }
}
