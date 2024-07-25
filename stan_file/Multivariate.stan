data { 
  int<lower=1> K; // number of respond variables
  //int<lower=1> J; // number of explained variables
  int<lower=0> N; // number of observations
  //vector[J] x[N];
  vector[N] x; // intercept term
  vector[K] y[N]; 
}
parameters { 
  //matrix[K, J] beta; 
  vector[K] beta; // intercept coefficient
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;
} 
model {
  //vector[K] mu[N];
  //for (n in 1:N) 
  //  mu[n] = beta * x[n];
  vector[K] mu[N];
  for (n in 1:N){
    mu[n] = beta * x[n];
  }
  
  
  //to_vector(beta) ~ normal(0, 2);
  beta ~ normal(0, 2);
  L_Omega ~ lkj_corr_cholesky(1); 
  L_sigma ~ student_t(3, 0, 2);
  
  y ~ multi_normal_cholesky(mu, diag_pre_multiply(L_sigma, L_Omega));
}
generated quantities {
  matrix[K, K] Omega;
  matrix[K, K] Sigma;
  Omega = multiply_lower_tri_self_transpose(L_Omega);
  Sigma = quad_form_diag(Omega, L_sigma); 
  
  vector[K] muhat[N];
  for(n in 1:N){
    muhat[n] = beta * x[n];
  }

  vector[K] Yhat[N];
  Yhat = multi_normal_cholesky_rng(muhat, diag_pre_multiply(L_sigma, L_Omega));
}
