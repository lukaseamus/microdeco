data{
  int n;
  vector[n] Dry;
}

parameters{
  real log_mu;
  real<lower=0> theta;
}

model{
  // Priors
  log_mu ~ normal( log(1.6) , 0.3 );
  theta ~ exponential( 1 );
  
  // Model
  real mu = exp( log_mu );
  Dry ~ gamma( mu / theta , 1 / theta );
}