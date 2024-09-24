data {
  int n;
  vector[n] disp;
  vector[n] climVelo;
  vector[n] shift;
}
transformed data {
  matrix[n,2] min_matrix = append_col(climVelo, disp);
  vector[n] min_rate;
 
  for(i in 1:n) {
    min_rate[i] = min(min_matrix[i,]);
  }
}
parameters {
 real b2;
 real<lower=0> sigma;
}
model {
  b2 ~ normal(1, 1);
  sigma ~ exponential(1);
  
  shift ~ gamma((min_rate*b2)^2 / sigma^2, (min_rate*b2) / sigma^2);
}
generated quantities {
  vector[n] shift_pred; 
  
  for(i in 1:n){
     shift_pred[i] = gamma_rng((min_rate[i]*b2)^2 / sigma^2, (min_rate[i]*b2) / sigma^2);
  }
}
