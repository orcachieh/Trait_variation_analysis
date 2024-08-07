data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  int V; // number of covariates
  vector[N] LW;   //measured leaf width
  matrix[N, V] CO; // covariates matix
}
transformed data{
}

parameters {
  positive_ordered[K] alpha;  // ordered intercept to identify
  
  simplex[K] theta; // mixture proportion
  
  array[K] simplex[V] beta_raw;
  real beta_scale;
  
  real<lower=0> sigma1;  // dispersion parameter
  real<lower=0> sigma2;  // dispersion parameter
}
transformed parameters{
  array[K] vector[V] beta;
  for(i in 1:K){
    beta[i] = beta_scale * (beta_raw[i] - inv(V));
  }
}

model {
  // priors
  for(i in 1:K){
  target += normal_lpdf(alpha[i] | 0, 2);
  target += beta_lpdf(beta_raw[i] | 1, 1);
  }
  target += normal_lpdf(beta_scale | 0, 5);
  
  target += lognormal_lpdf(sigma1 | 0, 2)- 
    normal_lccdf(0 | 0, 1);
  target += lognormal_lpdf(sigma2 | 0, 2)- 
    normal_lccdf(0 | 0, 1);
  
  target += beta_lpdf(theta | 1, 1);

  
  // likelihood

  for(n in 1:N){
    target += log_sum_exp(log(theta[1]) + 
                              lognormal_lpdf(LW[n] | alpha[1] + CO[n, ]*  beta[1], sigma1), // + bsed[1] * Sed_merge[n]
                            log(theta[2]) +
                              lognormal_lpdf(LW[n] | alpha[2] + CO[n, ] * beta[2], sigma2)); //+ bsed[2] * Sed_merge[n]
  }
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
      //      LW_gen[n] = lognormal_rng(alpha[1] + rei[n] * brei1 + de[n] * bde1, sigma1);
      //    } else {
        //      LW_gen[n] = lognormal_rng(alpha[2] + rei[n] * brei2 + de[n] * bde2, sigma2);
        //    }
    //  }
  //}

