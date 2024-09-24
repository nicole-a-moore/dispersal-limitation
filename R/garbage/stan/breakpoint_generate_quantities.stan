//can change names of data
data {
  int n;
  vector[n] disp;
  vector[n] climVelo;
}
// transformed data {
//   matrix[n,2] min_matrix = append_col(climVelo, disp);
//   vector[n] min_rate;
//  
//   for(i in 1:n) {
//     min_rate[i] = min(min_matrix[i,]);
//   }
// }
//don't mess with this
parameters {
 real b2;
 real<lower=0> sigma;
}
//delete model block
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
