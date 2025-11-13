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
  vector[n] Proportion_mean;
  vector[n] Proportion_sd;
  array[n] int Treatment;
  int n_Treatment;
  array[n] int Tank;
  int n_Tank;
}

transformed data{
  // Convert sd to nu because this is easier on the sampler
  vector[n] Proportion_nu = 
  Proportion_mean .* ( 1 + Proportion_mean ) ./ Proportion_sd^2;
}

parameters{
  // Parameter describing true, unobserved proportion
  vector<lower=0>[n] p;
  
  // Parameters describing mean
  /// Global parameters
  real log_alpha_mu;
  real log_mu_mu;
  real log_tau_mu;
  
  real<lower=0> log_alpha_sigma_t;
  real<lower=0> log_mu_sigma_t;
  real<lower=0> log_tau_sigma_t;
  
  real<lower=0> log_alpha_sigma_ta;
  real<lower=0> log_mu_sigma_ta;
  real<lower=0> log_tau_sigma_ta;
  
  /// Treatment parameters
  vector[n_Treatment] log_alpha_z_t;
  vector[n_Treatment] log_mu_z_t;
  vector[n_Treatment] log_tau_z_t;
  
  /// Tank parameters
  vector[n_Tank] log_alpha_z_ta;
  vector[n_Tank] log_mu_z_ta;
  vector[n_Tank] log_tau_z_ta;
  
  // Parameters describing precision
  /// Global parameters
  real log_epsilon_mu;
  real log_lambda_mu;
  real log_theta_mu;

  real<lower=0> log_epsilon_sigma_t;
  real<lower=0> log_lambda_sigma_t;
  real<lower=0> log_theta_sigma_t;
  
  real<lower=0> log_epsilon_sigma_ta;
  real<lower=0> log_lambda_sigma_ta;
  real<lower=0> log_theta_sigma_ta;
  
  /// Treatment parameters
  vector[n_Treatment] log_epsilon_z_t;
  vector[n_Treatment] log_lambda_z_t;
  vector[n_Treatment] log_theta_z_t;
  
  /// Tank parameters
  vector[n_Tank] log_epsilon_z_ta;
  vector[n_Tank] log_lambda_z_ta;
  vector[n_Tank] log_theta_z_ta;
}

model{
  // Priors
  /// Likelihood mean
  //// Global parameters
  log_alpha_mu ~ normal( log(0.001) , 0.3 );
  log_mu_mu ~ normal( log(40) , 0.3 );
  log_tau_mu ~ normal( log(0.06) , 0.3 );
  
  log_alpha_sigma_t ~ normal( 0 , 0.3 ) T[0,]; // half-normal prior
  log_mu_sigma_t ~ normal( 0 , 0.3 ) T[0,];
  log_tau_sigma_t ~ normal( 0 , 0.3 ) T[0,];
  
  log_alpha_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  log_mu_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  log_tau_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  
  //// Treatment parameters
  log_alpha_z_t ~ normal( 0 , 1 );
  log_mu_z_t ~ normal( 0 , 1 );
  log_tau_z_t ~ normal( 0 , 1 );
  
  vector[n_Treatment] log_alpha_t = log_alpha_z_t * log_alpha_sigma_t + log_alpha_mu;
  vector[n_Treatment] log_mu_t = log_mu_z_t * log_mu_sigma_t + log_mu_mu;
  vector[n_Treatment] log_tau_t = log_tau_z_t * log_tau_sigma_t + log_tau_mu;
  
  //// Tank parameters
  log_alpha_z_ta ~ normal( 0 , 1 );
  log_mu_z_ta ~ normal( 0 , 1 );
  log_tau_z_ta ~ normal( 0 , 1 );
  
  vector[n_Tank] log_alpha_ta = log_alpha_z_ta * log_alpha_sigma_ta + 0;
  vector[n_Tank] log_mu_ta = log_mu_z_ta * log_mu_sigma_ta + 0;
  vector[n_Tank] log_tau_ta = log_tau_z_ta * log_tau_sigma_ta + 0;
  
  /// Likelihood precision
  //// Global parameters
  log_epsilon_mu ~ normal( log(4e4) , 0.3 );
  log_lambda_mu ~ normal( log(0.1) , 0.3 );
  log_theta_mu ~ normal( log(500) , 0.3 );
  
  log_epsilon_sigma_t ~ normal( 0 , 0.3 ) T[0,];
  log_lambda_sigma_t ~ normal( 0 , 0.3 ) T[0,];
  log_theta_sigma_t ~ normal( 0 , 0.3 ) T[0,];
  
  log_epsilon_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  log_lambda_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  log_theta_sigma_ta ~ normal( 0 , 0.3 ) T[0,];
  
  //// Treatment parameters
  log_epsilon_z_t ~ normal( 0 , 1 );
  log_lambda_z_t ~ normal( 0 , 1 );
  log_theta_z_t ~ normal( 0 , 1 );
  
  vector[n_Treatment] log_epsilon_t = log_epsilon_z_t * log_epsilon_sigma_t + log_epsilon_mu;
  vector[n_Treatment] log_lambda_t = log_lambda_z_t * log_lambda_sigma_t + log_lambda_mu;
  vector[n_Treatment] log_theta_t = log_theta_z_t * log_theta_sigma_t + log_theta_mu;
  
  //// Tank parameters
  log_epsilon_z_ta ~ normal( 0 , 1 );
  log_lambda_z_ta ~ normal( 0 , 1 );
  log_theta_z_ta ~ normal( 0 , 1 );
  
  vector[n_Tank] log_epsilon_ta = log_epsilon_z_ta * log_epsilon_sigma_ta + 0;
  vector[n_Tank] log_lambda_ta = log_lambda_z_ta * log_lambda_sigma_ta + 0;
  vector[n_Tank] log_theta_ta = log_theta_z_ta * log_theta_sigma_ta + 0;
  
  // Model
  /// Likelihood mean
  //// Parameters
  vector[n] alpha = exp( log_alpha_t[Treatment] + log_alpha_ta[Tank] );
  vector[n] mu = exp( log_mu_t[Treatment] + log_mu_ta[Tank] );
  vector[n] tau = exp( log_tau_t[Treatment] + log_tau_ta[Tank] );
  
  //// Function
  vector[n] p_mu = exp(
      Day .* alpha - ( alpha + tau ) .* 
      mu ./ 5 .* (
        log1p_exp( 5 ./ mu .* ( Day - mu ) ) - log1p_exp( -5 )
      )
    );
  
  /// Likelihood precision
  //// Parameters
  vector[n] epsilon = exp( log_epsilon_t[Treatment] + log_epsilon_ta[Tank] );
  vector[n] lambda = exp( log_lambda_t[Treatment] + log_lambda_ta[Tank] );
  vector[n] theta = exp( log_theta_t[Treatment] + log_theta_ta[Tank] );
  
  //// Function
  vector[n] nu = theta + ( epsilon - theta ) .* exp(-lambda .* Day);
  
  // Beta prime likelihood
  for ( i in 1:n ) { // loop because betap isn't vectorised
    p[i] ~ betap( p_mu[i] * ( 1 + nu[i] ) , 2 + nu[i] );
  }
  
  // Beta prime measurement error model
  for ( i in 1:n ) {
    Proportion_mean[i] ~ betap(
      p[i] * ( 1 + Proportion_nu[i] ),
      2 + Proportion_nu[i]
    );
  }
  // for ( i in 1:n ) {
  //   Proportion_mean[i] ~ betap(
  //     p[i] * ( 1 + p[i] * ( 1 + p[i] ) / Proportion_sd[i]^2 ),
  //     2 + p[i] * ( 1 + p[i] ) / Proportion_sd[i]^2
  //   );
  // }
}

generated quantities{
  vector[n_Treatment] log_alpha_t = log_alpha_z_t * log_alpha_sigma_t + log_alpha_mu;
  vector[n_Treatment] log_mu_t = log_mu_z_t * log_mu_sigma_t + log_mu_mu;
  vector[n_Treatment] log_tau_t = log_tau_z_t * log_tau_sigma_t + log_tau_mu;
  
  vector[n_Tank] log_alpha_ta = log_alpha_z_ta * log_alpha_sigma_ta + 0;
  vector[n_Tank] log_mu_ta = log_mu_z_ta * log_mu_sigma_ta + 0;
  vector[n_Tank] log_tau_ta = log_tau_z_ta * log_tau_sigma_ta + 0;
  
  vector[n_Treatment] log_epsilon_t = log_epsilon_z_t * log_epsilon_sigma_t + log_epsilon_mu;
  vector[n_Treatment] log_lambda_t = log_lambda_z_t * log_lambda_sigma_t + log_lambda_mu;
  vector[n_Treatment] log_theta_t = log_theta_z_t * log_theta_sigma_t + log_theta_mu;
  
  vector[n_Tank] log_epsilon_ta = log_epsilon_z_ta * log_epsilon_sigma_ta + 0;
  vector[n_Tank] log_lambda_ta = log_lambda_z_ta * log_lambda_sigma_ta + 0;
  vector[n_Tank] log_theta_ta = log_theta_z_ta * log_theta_sigma_ta + 0;
}