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
 real<lower=0> sigma;
}
model {
  // add sampling error around dispersal and climate velocity observations 
  sigma ~ exponential(1);
  
  for(i in 1:n){
    if (disp[i] < climVelo[i]) {
      shift[i] ~ gamma((disp[i])^2 / sigma^2, (disp[i]) / sigma^2);
    } else {
      shift[i] ~ gamma((climVelo[i])^2 / sigma^2, (climVelo[i]) / sigma^2);
    }
  }
  
}
generated quantities {
  vector[n] shift_pred; // posterior predictive distribution
  // model's predictions (one per sample) predict what shifts should be given model
  // compare these to real data (posterior predictive check - step 1 in model validation)
  for(i in 1:n){
    if (disp[i] < climVelo[i]) {
      shift_pred[i] = gamma_rng((disp[i])^2 / sigma^2, (disp[i]) / sigma^2);
    } else {
      shift_pred[i] = gamma_rng((climVelo[i])^2 / sigma^2, (climVelo[i]) / sigma^2);
    }
  }
}
