data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  vector[N] LW;   //measured leaf width
  vector[N] rei;
}

parameters {
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
  target += normal_lpdf(brei1 | 0, 1);
  target += normal_lpdf(brei2 | 0, 1);
  
  target += normal_lpdf(alpha | 1, 0.5);
  target += normal_lpdf(gamma | 1, 0.5) - 
    normal_lcdf(alpha | 1, 0.5);
  
  target += lognormal_lpdf(sigma1 | 0, 1) - 
    normal_lccdf(0 | 0, 1);
  target += lognormal_lpdf(sigma2 | 0, 1) -
    normal_lccdf(0 | 0, 1);
  
  target += beta_lpdf(theta1 | 1, 1);
  target += beta_lpdf(theta2 | 1, 1);
  
  // likelihood
  for(n in 1:N)
    target += log_sum_exp(log(theta1) +
                            lognormal_lpdf(LW[n] | alpha + rei[n] * brei1, sigma1),
                          log(theta2) +
                            lognormal_lpdf(LW[n] | gamma + rei[n] * brei2, sigma2));
}

generated quantities{
  // real mu for each group
  //real mu1 = alpha + dot_product(rei, brei1);
  //real mu2 = gamma + dot_product(rei, brei2);
  
  // we will generate a component identifier "comp"
  // then draw "ytilde" using correct component mu and sigma
  int<lower=1,upper=K> comp[N];
  vector[N] LW_gen;
  
  //can draw N new observations or as many as you want
  for (n in 1:N) {
    comp[n] = categorical_rng(theta);
    if(comp[n] == 1){
    LW_gen[n] = lognormal_rng(alpha + rei[n] * brei1, sigma1);
    } else {
    LW_gen[n] = lognormal_rng(gamma + rei[n] * brei2, sigma2);
    }
  }
}
