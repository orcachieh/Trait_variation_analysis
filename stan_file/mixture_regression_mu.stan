data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  array[N] real LW;   //measured leaf width
}

parameters {
  simplex[K] theta; // mixture proportion
  vector<lower=0>[K] sigma;  // dispersion parameter
  ordered[K] mu;
  }

model {
  // priors
  
  vector[K] log_theta = log(theta);  // cache log calculation
  sigma ~ lognormal(0, 2);
  mu ~ normal(0, 10);

  for (n in 1:N) {
    vector[K] lps = log_theta;
    for (k in 1:K)
    lps[k] += normal_lpdf(LW[n] | mu[k], sigma[k]);
    target += log_sum_exp(lps);
  }
}
