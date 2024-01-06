data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  vector[N] LW;   //measured leaf width
  vector[N] de;
}

parameters {
  ordered[K] b;
  
  simplex[K] theta; // mixture proportion
  vector<lower=0>[K] sigma;  // dispersion parameter
  
}
transformed parameters {
  ordered[2] mu;
  for(k in 1:K){
    mu[k] = mean(b[k] * de);
  }
}

model {
  // priors
  b ~ normal(1 , 0.5); 
  
  vector[K] log_theta = log(theta);  // cache log calculation
  sigma ~ lognormal(0, 2);

  for (n in 1:N) {
    vector[K] lps = log_theta;
    for(k in 1:K)
    lps[k] += normal_lpdf(LW[n] | mu[k], sigma[k]);
    target += log_sum_exp(lps);
  }
}

//generated quantities {
  // actual population-level intercept
  //real b_mu1_Intercept = Intercept_mu1 - dot_product(X1, beta1);
  // actual population-level intercept
  //real b_mu2_Intercept = Intercept_mu2 - dot_product(X2, beta2);
  //}

