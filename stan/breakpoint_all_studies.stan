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
 real b2;
 real<lower=0> sigma;
}
model {
  b2 ~ normal(1, 1);
  sigma ~ exponential(1);
  
  for(i in 1:n){
    if (disp[i] < climVelo[i]) {
      shift[i] ~ gamma((disp[i]*b2)^2 / sigma^2, (disp[i]*b2) / sigma^2);
    } else {
      shift[i] ~ gamma((climVelo[i]*b2)^2 / sigma^2, (climVelo[i]*b2) / sigma^2);
    }
  }
  
}
generated quantities {
  vector[n] shift_pred; // posterior predictive distribution
  // model's predictions (one per sample) predict what shifts should be given model
  // compare these to real data (posterior predictive check - step 1 in model validation)
  for(i in 1:n){
    if (disp[i] < climVelo[i]) {
      shift_pred[i] = gamma_rng((disp[i]*b2)^2 / sigma^2, (disp[i]*b2) / sigma^2);
    } else {
      shift_pred[i] = gamma_rng((climVelo[i]*b2)^2 / sigma^2, (climVelo[i]*b2) / sigma^2);
    }
  }
}
