data{
  int n;
  vector[n] Day;
  vector[n] GP_P;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int Tank;
  int n_Tank;
}

parameters{
  // Likelihood mean
  /// Global parameters
  real log_mu_mu;
  real<lower=0> log_mu_sigma_t;
  real<lower=0> log_mu_sigma_ta;
  
  /// Treatment parameters
  vector[n_Treatment] log_mu_t;
  
  /// Tank parameters
  vector[n_Tank] log_mu_ta;

  // Likelihood precision
  real<lower=0> nu;
}

model{
  // Priors
  /// Likelihood mean
  //// Global parameters
  log_mu_mu ~ normal( log(50) , 0.3 );
  log_mu_sigma_t ~ normal( 0 , 0.3 ) T[0,];
  log_mu_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  
  //// Treatment parameters
  log_mu_t ~ normal( log_mu_mu , log_mu_sigma_t );
  
  //// Tank parameters
  log_mu_ta ~ normal( 0 , log_mu_sigma_ta );
  
  /// Likelihood precision
  nu ~ gamma( square(30) / square(20) , 30 / square(20) );

  // Model
  /// Likelihood mean
  //// Parameter
  vector[n] mu = exp( log_mu_t[Treatment] + log_mu_ta[Tank] );
  
  //// Function
  vector[n] P_mu = inv_logit( -5 / mu .* ( Day - mu ) );
  
  /// Beta likelihood
  GP_P ~ beta( P_mu * nu , (1 - P_mu) * nu );
}