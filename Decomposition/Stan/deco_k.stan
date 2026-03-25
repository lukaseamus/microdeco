data{
  int n;
  vector[n] Day;
  vector[n] Ratio_mean;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int Tank;
  int n_Tank;
}

parameters{
  real log_k_mu;
  real<lower=0> log_k_sigma_t;
  real<lower=0> log_k_sigma_ta;
  real<lower=0> sigma;

  vector[n_Treatment] log_k_z_t;
  vector[n_Tank] log_k_z_ta;
}

transformed parameters{
  vector[n_Treatment] log_k_t = log_k_z_t * log_k_sigma_t + log_k_mu;
  vector[n_Tank] log_k_ta = log_k_z_ta * log_k_sigma_ta + 0;
}

model{
  // Priors
  log_k_mu ~ normal( log(0.06) , 1 );
  log_k_sigma_t ~ normal( 0 , 1 ) T[0,];
  log_k_sigma_ta ~ normal( 0 , 1 ) T[0,];
  sigma ~ exponential( 1 );
  
  log_k_z_t ~ normal( 0 , 1 );
  log_k_z_ta ~ normal( 0 , 1 );
  
  // Model
  vector[n] k = exp( log_k_t[Treatment] + log_k_ta[Tank] );
  vector[n] r_mu = exp( -k .* Day );
  
  // Normal likelihood
  Ratio_mean ~ normal( r_mu , sigma );
}