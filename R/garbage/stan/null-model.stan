data {
  int n;
  vector[n] disp;
  vector[n] climVelo;
  vector[n] shift;
}
// transformed data {
//   real M = max(x);
// }
parameters {
 real<lower=0> mu;
 real<lower=0> sigma;
}
model {
  mu ~ normal(1, 1);
  sigma ~ exponential(1);
  // vector[n] mu = alpha0 + alpha
  
  shift ~ gamma(mu^2/sigma^2, mu/sigma^2);
}
generated quantities {
  vector[n] shift_pred; 
  for(i in 1:n){
      shift_pred[i] = gamma_rng((mu)^2 / sigma^2, (mu) / sigma^2);
  }
  
  vector[n] log_lik;
  for (j in 1:n) {
    log_lik[j] = gamma_lpdf(shift_pred[j] | (mu)^2 / sigma^2, (mu) / sigma^2);
    //log probability density function
    // what is the probability of each shift_pred observation given a gamma distribution with parameters fitted by the model 
  }
  
}
