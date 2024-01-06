data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  int P;          // intercept + number of predictors
  vector[N] LW;   //measured leaf width
  matrix[N, P] X1; // Predictor matrix
  matrix[N, P] X2; // Predictor matrix
}

parameters {
  vector[P] beta1; // coefficients set1
  vector[P] beta2; // coefficients set2
  
  simplex[K] theta; // mixture proportion
  vector<lower=0>[K] sigma;  // dispersion parameter
  ordered[2] ordered_Intercept;  // to identify mixtures
}
transformed parameters {
  // identify mixtures via ordering of the intercepts
  real Intercept_mu1 = ordered_Intercept[1];
  // identify mixtures via ordering of the intercepts
  real Intercept_mu2 = ordered_Intercept[2];
}

model {
  // priors
  beta1 ~ normal(1 , 0.5);     
  beta2 ~ normal(1, 0.5);      

  vector[K] log_theta = log(theta);  // cache log calculation
  sigma ~ lognormal(0, 2);
  
  vector[N] mu1 = rep_vector(0.0, N);
  vector[N] mu2 = rep_vector(0.0, N);
  mu1 += Intercept_mu1 + X1 * beta1;
  mu2 += Intercept_mu2 + X2 * beta2;
  
  for (n in 1:N) {
    vector[K] lps = log_theta;
    lps[1] = normal_lpdf(LW[n] | mu1[n], sigma[1]);
    lps[2] = normal_lpdf(LW[n] | mu2[n], sigma[2]);
    target += log_sum_exp(lps);
    }
}

//generated quantities {
  // actual population-level intercept
  //real b_mu1_Intercept = Intercept_mu1 - dot_product(X1, beta1);
  // actual population-level intercept
  //real b_mu2_Intercept = Intercept_mu2 - dot_product(X2, beta2);
//}

  