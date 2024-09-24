data {
  int n;
  int num_study; // number of studies 
  vector[n] disp;
  vector[n] climVelo;
  vector[n] shift;
  array[n] int<lower=1,upper=num_study> study_id; //unique study ids (1 - num_study)
}
// transformed data {
//   real M = max(x);
// }
parameters {
 real log_alpha;
 real<lower=0> sigma;
 vector[num_study] log_study_diffs; // vector of length num_study specifying the log of the study-level differences in the mean shift 
 real<lower=0> sigma_log_study_diffs; // the standard deviation of log study level differences in the mean shift 
}
model {
  log_alpha ~ normal(0, 1);
  sigma ~ exponential(1);
  sigma_log_study_diffs ~ exponential(1);
  log_study_diffs ~ normal(0, sigma_log_study_diffs); // hiercharchical model - uses data to inform priors

  // calculate study-level means
  vector[n] log_mu = log_alpha + log_study_diffs[study_id];
  // unlog study-level means
  vector[n] mu = exp(log_mu);
  // calculate shift
  shift ~ gamma(mu^2/sigma^2, mu/sigma^2);
}
generated quantities {
  vector[n] shift_pred;

  // calculate study-level means
  vector[n] log_mu = log_alpha + log_study_diffs[study_id];
  // unlog study-level means
  vector[n] mu = exp(log_mu);
  
  for(i in 1:n){
      shift_pred[i] = gamma_rng((mu[i])^2 / sigma^2, (mu[i]) / sigma^2);
  }
  
  vector[n] log_lik;
  for (j in 1:n) {
    log_lik[j] = gamma_lpdf(shift_pred[j] | (mu[j])^2 / sigma^2, (mu[j]) / sigma^2);
    //log probability density function
    // what is the probability of each shift_pred observation given a gamma distribution with parameters fitted by the model 
  }
}
