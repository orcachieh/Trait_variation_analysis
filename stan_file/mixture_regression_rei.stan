data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  vector[N] LW;   //measured leaf width
  vector[N] rei;
}

parameters {
  //real brei1;
  real brei2;
  
  positive_ordered[K] alpha;
  
  simplex[K] theta; // mixture proportion
  real<lower=0> sigma1;  // dispersion parameter
  real<lower=0> sigma2;  // dispersion parameter
  
}
transformed parameters {
}

model {
  // priors
  //target += normal_lpdf(brei1 | 0, 1);
  target += normal_lpdf(brei2 | 0, 1);
  
  target += normal_lpdf(alpha | 1, 0.5);

  target += lognormal_lpdf(sigma1 | 0, 1) - 
    normal_lccdf(0 | 0, 1);
  target += lognormal_lpdf(sigma2 | 0, 1) -
    normal_lccdf(0 | 0, 1);
  
  target += beta_lpdf(theta | 1, 1);
  
  // likelihood
  for(n in 1:N)
    target += log_sum_exp(log(theta[1]) +
                            normal_lpdf(LW[n] | alpha[1], sigma1),
                          log(theta[2]) +
                            normal_lpdf(LW[n] | alpha[2] + rei[n] * brei2, sigma2));
}

//generated quantities{
  // real mu for each group
  //real mu1 = alpha + dot_product(rei, brei1);
  //real mu2 = gamma + dot_product(rei, brei2);
  
  // we will generate a component identifier "comp"
  // then draw "ytilde" using correct component mu and sigma
//  int<lower=1,upper=K> comp[N];
//  vector[N] LW_gen;
  
  //can draw N new observations or as many as you want
//  for (n in 1:N) {
//    comp[n] = categorical_rng(theta);
//    if(comp[n] == 1){
//    LW_gen[n] = lognormal_rng(alpha + rei[n] * brei1, sigma1);
//    } else {
//    LW_gen[n] = lognormal_rng(gamma + rei[n] * brei2, sigma2);
//    }
//  }
//}
