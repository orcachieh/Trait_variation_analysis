data {
  int<lower=1> K;  //number of mixture component (clusters) 
  int<lower=1> N;  //      number of observations
  vector[N] LW;   //measured leaf width
  vector[N] D;    //depth
  vector[N] REI;  //REI
  array[N] int<lower=1> AIR;  //Air exposure (Ordered categories)
  vector[4] alpha;   //air exposure prior
  vector[N] GA;   //Gravel
  vector[N] SA;   //Sand
  vector[N] SI;   //Silt
  vector[N] CL;   //Clay
}

parameters {
  vector[K] d;           //for Depth
  vector[K] rei;         //for REI
  vector[K] air;         //for Air exposure
  vector[K] ga;          //for Gravel
  vector[K] sa;          //for Sand
  vector[K] si;          //for Silt
  vector[K] cl;          //for Clay
  
  simplex[4] delta;  // a simplex to deal with ordered expained variable
  
  ordered[K] ordered_Intercept;  // to identify mixtures
  
  simplex[K] theta; // mixture proportion
  real<lower=0> sigma1;  // dispersion parameter
  real<lower=0> sigma2;  // dispersion parameter
}
transformed parameters{
  // identify mixtures via ordering of the intercepts
  real Intercept_mu1 = ordered_Intercept[1];
  // identify mixtures via ordering of the intercepts
  real Intercept_mu2 = ordered_Intercept[2];
  
  real d1;
  d1 = d[1];
  real d2;
  d2 = d[2];
  real rei1;
  rei1 = rei[1];
  real rei2;
  rei2 = rei[2];
  real air1;
  air1 = air[1];
  real air2;
  air2 = air[2];
  real ga1;
  ga1 = ga[1];
  real ga2;
  ga2 = ga[2];
  real sa1;
  sa1 = sa[1];
  real sa2;
  sa2 = sa[2];
  real si1;
  si1 = si[1];
  real si2;
  si2 = si[2];
  real cl1;
  cl1 = cl[1];
  real cl2;
  cl2 = cl[2];
  
  real<lower=0,upper=1> theta1;
  real<lower=0,upper=1> theta2;
  theta1 = theta[1];
  theta2 = theta[2];
}
model {
  // priors
  d1 ~ normal(1 , 0.5);     
  rei1 ~ normal(1, 0.5);      
  air1 ~ normal(1, 0.5);
  ga1 ~ normal(1, 0.5);
  sa1 ~ normal(1, 0.5);
  si1 ~ normal(1, 0.5);
  cl1 ~ normal(1, 0.5);
  d2 ~ normal(1 , 0.5);     
  rei2 ~ normal(1, 0.5);    
  air2 ~ normal(1, 0.5);
  ga2 ~ normal(1, 0.5);
  sa2 ~ normal(1, 0.5);
  si2 ~ normal(1, 0.5);
  cl2 ~ normal(1, 0.5);
  
  vector[5] delta_j;
  delta ~ dirichlet(alpha);
  delta_j = append_row(0, delta);

  sigma1 ~ lognormal(0, 2);
  sigma2 ~ lognormal(0, 2);
  
  theta1 ~ lognormal(0, 0.1);
  theta2 ~ lognormal(0, 0.1);
  
  vector[N] mu1 = rep_vector(0.0, N);
  vector[N] mu2 = rep_vector(0.0, N);
  
  for(n in 1:N){
      mu1[n] = Intercept_mu1 + d1 * D[n] + rei1 * REI[n] + air1 * sum(delta_j[1:AIR[n]]) + ga1 * GA[n] + sa1 * SA[n] + si1 * SI[n] + cl1 * CL[n];
      mu2[n] = Intercept_mu2 + d2 * D[n] + rei2 * REI[n] + air2 * sum(delta_j[1:AIR[n]]) + ga2 * GA[n] + sa2 * SA[n] + si2 * SI[n] + cl2 * CL[n];
  }
  
//  for(n in 1:N){
//    mu1[n] = Intercept_mu1 + d1 * D[n] + rei1 * REI[n];
//    mu2[n] = Intercept_mu2 + d2 * D[n] + rei2 * REI[n];
//  }

//   mu[K] = d[K] * D + rei[K] * REI + air[K] * sum(air9[1:AIR]) + ga[K] * GA + sa[K] * SA + si[K] * SI + cl[K] * CL;
  
  
  for (n in 1:N) {
    real ps[2];
    ps[1] = log(theta1) + normal_lpdf(LW[n] | mu1[n], sigma1);
    ps[2] = log(theta2) + normal_lpdf(LW[n] | mu2[n], sigma2);
    target += log_sum_exp(ps);
  }
}

//generated quantities {
  // we will generate a component identifier "comp"
  // then draw "ytilde" using correct component mu and sigma
//  int<lower=1,upper=K> comp[N];
//  vector[N] ytilde;
  
  // can draw N new observations or as many as you want
//  for (n in 1:N) {
//    comp[n] = categorical_rng(theta);
//    if(comp[n] == 1){
//    ytilde[n] = normal_rng(mu[comp[n]], sigma[comp[n]]);
//    } else {
//    ytilde[n] = normal_rng(mu[comp[n]], sigma[comp[n]]);
//    }

//  matrix[K,N] membership_probs ;
//  {
//    vector[K] log_theta = log(theta);  // cache log calculation
//    for (n in 1:N) {
//      vector[K] lps = log_theta;
//      for (k in 1:K)
//        lps[k] += normal_lpdf(y[n] | mu[k], sigma[k]);
//      membership_probs[,n] = exp(lps) ;
//    }
  }
//}
