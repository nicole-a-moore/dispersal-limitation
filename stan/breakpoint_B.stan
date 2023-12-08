data {
  int n;
  vector[n] disp;
  real climVelo;
  vector[n] shift;
}
// transformed data {
//   real M = max(x);
// }
parameters {
 real b2;
 real<lower=0> sigma;
}
model {
  b2 ~ normal(1, 1);
  sigma ~ exponential(1);
  
  for(i in 1:n){
    if (disp[i] < climVelo) {
      shift[i] ~ normal(disp[i]*b2, sigma);
    } else {
      shift[i] ~ normal(climVelo*b2, sigma);
    }
  }
  
}
generated quantities {
  vector[n] shift_pred; // posterior predictive distribution
  // model's predictions (one per sample) predict what shifts should be given model
  // compare these to real data (posterior predictive check - step 1 in model validation)
  for(i in 1:n){
    if (disp[i] < climVelo) {
      shift_pred[i] = normal_rng(disp[i]*b2, sigma);
    } else {
      shift_pred[i] = normal_rng(climVelo*b2, sigma);
    }
  }
}
