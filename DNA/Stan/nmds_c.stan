data{
  int n;
  int k;
  array[n] vector [k] y; 
  array[n] int Treatment;
  int n_Treatment;
}

parameters{
  // Likelihood mean
  vector[k] mu_mu; // Global mean
  vector<lower=0>[k] mu_sigma; // Standard deviation across treatments
  array[n_Treatment] vector[k] mu; // Treatment means

  // Likelihood covariance matrix
  corr_matrix[k] Omega; // Correlation matrix
  vector<lower=0>[k] sigma; // Standard deviations
}

transformed parameters{
  // Calculate likelihood covariance matrix
  matrix[k, k] Sigma = quad_form_diag( Omega , sigma );
}

model{
  // Priors
  /// Likelihood mean
  mu_mu ~ std_normal(); // Standard normal because y is standardised
  mu_sigma ~ std_normal() T[0,]; // Half-normal
  for (i in 1:n_Treatment) {
    mu[i] ~ normal( mu_mu , mu_sigma ); // Partial pooling
  }
  
  /// Likelihood covariance matrix
  Omega ~ lkj_corr( 2 );
  sigma ~ exponential( 1 );
  
  // Multivariate normal likelihood
  y ~ multi_normal( mu[Treatment] , Sigma );
}