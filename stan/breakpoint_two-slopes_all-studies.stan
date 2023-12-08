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
 real b1;
 real b2;
 real<lower=0> sigma;
}
model {
  b1 ~ normal(1, 1);
  b2 ~ normal(0, 1);
  sigma ~ exponential(1);
  
  for(i in 1:n){
    if (disp[i] < climVelo[i]) {
      shift[i] ~ gamma((disp[i]*b1)^2 / sigma^2, (disp[i]*b1) / sigma^2);
    } else {
      shift[i] ~ gamma((climVelo[i]*b1*b2)^2 / sigma^2, (climVelo[i]*b1*b2) / sigma^2);
    }
  }
  
}
generated quantities {
  vector[n] shift_pred; // posterior predictive distribution
  // model's predictions (one per sample) predict what shifts should be given model
  // compare these to real data (posterior predictive check - step 1 in model validation)
  for(i in 1:n){
    if (disp[i] < climVelo[i]) {
      shift_pred[i] = gamma_rng((disp[i]*b1)^2 / sigma^2, (disp[i]*b1) / sigma^2);
    } else {
      shift_pred[i] = gamma_rng((climVelo[i]*b1*b2)^2 / sigma^2, (climVelo[i]*b1*b2) / sigma^2);
    }
  }
}
