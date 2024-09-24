data {
  int n;
  int num_study; // number of studies 
  int num_sp; // number of spp 
  vector[n] disp;
  vector[n] climVelo;
  vector[n] shift;
  array[n] int<lower=1,upper=num_study> study_id; //unique study ids (1 - num_study)
  array[n] int<lower=1,upper=num_sp> sp_id; //unique species ids (1 - num_sp)
}
// transformed data {
//   real M = max(x);
// }
parameters {
 real log_alpha;
 real<lower=0> sigma;
 
 vector[num_study] log_study_diffs; // vector of length num_study specifying the log of the study-level differences in the mean shift 
 real<lower=0> sigma_log_study_diffs; // the standard deviation of log study level differences in the mean shift 
 
  vector[num_sp] log_sp_diffs; // vector of length num_study specifying the log of the species-level differences in the mean shift 
 real<lower=0> sigma_log_sp_diffs; // the standard deviation of log spp level differences in the mean shift 
}
model {
  log_alpha ~ normal(0, 1);
  sigma ~ exponential(1);
  
  sigma_log_study_diffs ~ exponential(1);
  log_study_diffs ~ normal(0, sigma_log_study_diffs); // hiercharchical model - uses data to inform priors
  
  sigma_log_sp_diffs ~ exponential(1);
  log_sp_diffs ~ normal(0, sigma_log_sp_diffs); 

  // calculate means
  vector[n] log_b2 = log_alpha + log_study_diffs[study_id] + log_sp_diffs[sp_id];
  // unlog means
  vector[n] b2 = exp(log_b2);
  
  // calculate shift
  for(i in 1:n){
    if (disp[i] < climVelo[i]) {
      shift[i] ~ gamma((disp[i]*b2[i])^2 / sigma^2, (disp[i]*b2[i]) / sigma^2);
    } else {
      shift[i] ~ gamma((climVelo[i]*b2[i])^2 / sigma^2, (climVelo[i]*b2[i]) / sigma^2);
    }
  }
}
generated quantities {
  vector[n] shift_pred;

  // calculate study-level means
  vector[n] log_b2 = log_alpha + log_study_diffs[study_id] + log_sp_diffs[sp_id];
  // unlog study-level means
  vector[n] b2 = exp(log_b2);
  
  for(i in 1:n){
    if (disp[i] < climVelo[i]) {
      shift_pred[i] = gamma_rng((disp[i]*b2[i])^2 / sigma^2, (disp[i]*b2[i]) / sigma^2);
    } else {
      shift_pred[i] = gamma_rng((climVelo[i]*b2[i])^2 / sigma^2, (climVelo[i]*b2[i]) / sigma^2);
    }
  }
}
