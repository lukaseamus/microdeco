functions{
  // Beta prime log probability density function
  real betap_lpdf( real y , real alpha , real beta ) {
    return ( alpha - 1 ) * log( y )
    - ( alpha + beta ) * log1p( y ) -
    lbeta( alpha , beta );
  }
}

data{
  int n;
  vector[n] Day;
  vector[n] Ratio_mean;
  vector[n] Ratio_sd;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int Tank;
  int n_Tank;
}

transformed data{
  // Convert sd to nu because this is easier on the sampler
  vector[n] Ratio_nu = Ratio_mean .* ( 1 + Ratio_mean ) ./ Ratio_sd^2;
}

parameters{
  // Parameter describing true, unobserved proportion
  vector<lower=0>[n] r;
  
  // Parameters describing mean
  /// Global parameters
  real log_delta_mu; // delta = alpha + tau
  real log_mu_mu;
  real log_tau_mu;
  
  real<lower=0> log_delta_sigma_t;
  real<lower=0> log_mu_sigma_t;
  real<lower=0> log_tau_sigma_t;
  
  real<lower=0> log_delta_sigma_ta;
  real<lower=0> log_mu_sigma_ta;
  real<lower=0> log_tau_sigma_ta;
  
  /// Treatment parameters
  vector[n_Treatment] log_delta_t;
  vector[n_Treatment] log_mu_t;
  vector[n_Treatment] log_tau_t;
  
  /// Tank parameters
  vector[n_Tank] log_delta_ta;
  vector[n_Tank] log_mu_ta;
  vector[n_Tank] log_tau_ta;
  
  // Parameters describing precision
  real<lower=0> epsilon;
  real<lower=0> lambda;
  real<lower=0> theta;
}

model{
  // Priors
  /// Likelihood mean
  //// Global parameters
  log_delta_mu ~ normal( log(0.05) , 0.3 );
  log_mu_mu ~ normal( log(50) , 0.3 );
  log_tau_mu ~ normal( log(0.06) , 0.3 );
  
  log_delta_sigma_t ~ normal( 0 , 0.3 ) T[0,]; // half-normal prior
  log_mu_sigma_t ~ normal( 0 , 0.3 ) T[0,];
  log_tau_sigma_t ~ normal( 0 , 0.3 ) T[0,];
  
  log_delta_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  log_mu_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  log_tau_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  
  //// Treatment parameters
  log_delta_t ~ normal( log_delta_mu , log_delta_sigma_t );
  log_mu_t ~ normal( log_mu_mu , log_mu_sigma_t );
  log_tau_t ~ normal( log_tau_mu , log_tau_sigma_t );
  
  //// Tank parameters
  log_delta_ta ~ normal( 0 , log_delta_sigma_ta );
  log_mu_ta ~ normal( 0 , log_mu_sigma_ta );
  log_tau_ta ~ normal( 0 , log_tau_sigma_ta );
  
  /// Likelihood precision
  epsilon ~ gamma( square(4e4) / square(2e4) , 4e4 / square(2e4) );
  lambda ~ exponential( 1 );
  theta ~ gamma( square(500) / square(250) , 500 / square(250) );

  // Model
  /// Likelihood mean
  //// Parameters
  vector[n] delta = exp( log_delta_t[Treatment] + log_delta_ta[Tank] );
  vector[n] mu = exp( log_mu_t[Treatment] + log_mu_ta[Tank] );
  vector[n] tau = exp( log_tau_t[Treatment] + log_tau_ta[Tank] );
  vector[n] alpha = delta - tau;
  
  //// Function
  vector[n] r_mu = exp(
    Day .* alpha - ( alpha + tau ) .* mu ./ 5 .* (
      log1p_exp( 5 ./ mu .* ( Day - mu ) ) - log1p_exp( -5 )
    )
  );
  
  /// Likelihood precision
  // This parameterisation is better for large variance (partial pooling on variance)
  // vector[n] nu = theta + ( epsilon - theta ) .* exp( -lambda .* Day );
  
  // This parameterisation is typically better
  vector[n] nu = theta + exp( log( epsilon - theta ) - lambda .* Day );
  
  // Beta prime likelihood
  for ( i in 1:n ) r[i] ~ betap( r_mu[i] * ( 1 + nu[i] ) , 2 + nu[i] );
  
  // Beta prime measurement error model
  for ( i in 1:n ) Ratio_mean[i] ~ betap( 
    r[i] * ( 1 + Ratio_nu[i] ),
    2 + Ratio_nu[i]
  );
  
  // Alternative mean-sd parameterisation
  // for ( i in 1:n ) {
  //   Ratio_mean[i] ~ betap(
  //     r[i] * ( 1 + r[i] * ( 1 + r[i] ) / Ratio_sd[i]^2 ),
  //     2 + r[i] * ( 1 + r[i] ) / Ration_sd[i]^2
  //   );
  // }
}