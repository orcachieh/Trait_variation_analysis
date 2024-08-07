functions{
  
  vector merge_missing( array[] int miss_indexes , vector x_obs , vector x_miss ) {
    int N = dims(x_obs)[1];
    int N_miss = dims(x_miss)[1];
    vector[N] merged;
    merged = x_obs;
    for ( i in 1:N_miss )
      merged[ miss_indexes[i] ] = x_miss[i];
    return merged;
  }
}
data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  vector[N] LW;   //measured leaf width
  vector[N] sed;   //  mean sediment size to inform regression in each mixture component. contains 10 missing value
  array[10] int Sed_missidx;
}

parameters {
  positive_ordered[K] alpha;  // ordered intercept to identify
  array[K] real bsed; // one coefficient for each component
  
  simplex[K] theta; // mixture proportion
  real<lower=0> sigma1;  // dispersion parameter
  real<lower=0> sigma2;  // dispersion parameter
  
  real nu;
  real sigma_S;
  
  vector[10] Sed_impute;
}
transformed parameters{
}

model {
  // priors
  for(i in 1:K){
    target += normal_lpdf(bsed[i] | 0, 0.5);
    target += normal_lpdf(alpha[i] | 0, 1);
  }
  
  target += normal_lpdf(sigma1 | 0, 0.5);
  target += normal_lpdf(sigma2 | 0, 0.5);
  
  target += beta_lpdf(theta[1] | 2, 0.7);
  target += beta_lpdf(theta[2] | 0.7, 2);
  
  
  // likelihood
  vector[N] Sed_merge;
  Sed_merge = merge_missing(Sed_missidx, to_vector(sed), Sed_impute);
  Sed_merge ~ normal(nu, sigma_S);
  
  for(n in 1:N){
    target += log_sum_exp(log(theta[1]) + 
                              normal_lpdf(LW[n] | alpha[1] + bsed[1] * Sed_merge[n], sigma1), 
                            log(theta[2]) +
                              normal_lpdf(LW[n] | alpha[2] + bsed[2] * Sed_merge[n] , sigma2));
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

