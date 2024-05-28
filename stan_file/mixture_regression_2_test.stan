data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  vector[N] LW;   //measured leaf width
  vector[N] de;
  vector[N] rei;
}

parameters {
  real bde1;
  real bde2;
  real brei1;
  real brei2;
  
  real alpha;
  real<upper = alpha> gamma;
  
  simplex[K] theta; // mixture proportion
  real<lower=0> sigma1;  // dispersion parameter
  real<lower=0> sigma2;  // dispersion parameter
  
}
transformed parameters {
  real<lower = 0, upper = 1> theta1;
  real<lower = 0, upper = 1> theta2;
  theta1 = theta[1];
  theta2 = theta[2];
}

model {
  // priors
  target += normal_lpdf(bde1 | 0, 1);
  target += normal_lpdf(bde2 | 0, 1);
  target += normal_lpdf(brei1 | 0, 1);
  target += normal_lpdf(brei2 | 0, 1);
  
  target += normal_lpdf(alpha | 6, 1);
  target += normal_lpdf(gamma | 6, 1) - 
   normal_lcdf(alpha | 6, 1);
  
  target += lognormal_lpdf(sigma1 | 0, 1) - 
   normal_lccdf(0 | 0, 1);
  target += lognormal_lpdf(sigma2 | 0, 1) -
   normal_lccdf(0 | 0, 1);
  
  target += beta_lpdf(theta1 | 1, 1);
  target += beta_lpdf(theta2 | 1, 1);
  
  // likelihood
  for(n in 1:N)
    target += log_sum_exp(log(theta1) +
                          normal_lpdf(LW[n] | alpha + de[n] * bde1 + rei[n] * brei1, sigma1),
                          log(theta2) +
                          normal_lpdf(LW[n] | gamma + de[n] * bde2 + rei[n] * brei2, sigma2));
}

//generated quantities {
  // actual population-level intercept
  //real b_mu1_Intercept = Intercept_mu1 - dot_product(X1, beta1);
  // actual population-level intercept
  //real b_mu2_Intercept = Intercept_mu2 - dot_product(X2, beta2);
  //}

