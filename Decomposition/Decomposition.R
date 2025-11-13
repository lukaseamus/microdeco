##### Decomposition ####
# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)
deco <- here("Decomposition", "Decomposition.csv") %>% 
  read_csv() %>%
  mutate(Date = Date %>% dmy(),
         Day = Date[1] %--% Date / ddays(),
         Tank = Tank %>% factor(),
         Treatment = case_when(
           PAR == 0 ~ "Dark 15°C",
           Temperature == 15 ~ "Light 15°C",
           Temperature == 20 ~ "Light 20°C",
           Temperature == 25 ~ "Light 25°C"
         ) %>% fct_relevel("Dark 15°C"),
         Dry = if_else(Dry == 0, 1e-3, Dry)) %T>%
  print(n = 122)

# 1.2 Estimate initial mass ####
# Call model
require(cmdstanr)
initial_model <- here("Decomposition", "Stan", "initial.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

# Run model
require(tidybayes)
initial_samples <- initial_model$sample(
  data = deco %>%
    filter(Day == 0) %>%
    select(Dry) %>%
    compose_data(),
  chains = 8,
  parallel_chains = parallel::detectCores(),
  iter_warmup = 1e4,
  iter_sampling = 1e4
) %T>%
  print() # Rhat and ESS look great

# Check chains
require(bayesplot)
initial_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

# Prior-posterior comparison
source("functions.R")
prior_samples(
  model = initial_model,
  data = deco %>%
    filter(Day == 0) %>%
    select(Dry) %>%
    compose_data()
) %>%
  prior_posterior_draws(
    posterior_samples = initial_samples,
    parameters = c("log_mu", "theta"),
    format = "long"
  ) %>%
  prior_posterior_plot()

# Prediction
initial_posterior <- initial_samples %>%
  spread_draws(log_mu, theta) %>%
  mutate(mu = exp( log_mu ),
         Initial = rgamma( n() , mu / theta , 1 / theta )) %T>%
  print()

initial_posterior %>%
  summarise(Initial_mean = mean(Initial),
            Initial_sd = sd(Initial),
            n = n())

# 1.3 Calculations ####
deco %<>%
  filter(Day != 0) %>%
  cross_join(
    initial_posterior %>% 
      select(starts_with("."), Initial)
  ) %T>%
  print()

deco %<>%
  mutate(Proportion = Dry / Initial,
         k = ( log(Initial) - log(Dry) ) / Day)

deco_summary <- deco %>%
  group_by(Date, Tank, Temperature, PAR, Wet, Dry, d15N, d13C,
           N, C, C_N, Day, Treatment) %>%
  summarise(Initial_mean = mean(Initial),
            Initial_sd = sd(Initial),
            Proportion_mean = mean(Proportion),
            Proportion_sd = sd(Proportion),
            Proportion_lwr = qi(Proportion, .width = 0.90)[1],
            Proportion_upr = qi(Proportion, .width = 0.90)[2],
            k_mean = mean(k),
            k_sd = sd(k),
            k_lwr = qi(k, .width = 0.90)[1],
            k_upr = qi(k, .width = 0.90)[2],
            n = n()) %>%
  ungroup() %T>%
  print(n = 122)

# 2. Analyse data ####
# 2.1 Data exploration ####
deco_summary %>%
  ggplot() +
  geom_pointrange(aes(Day, Proportion_mean, 
                      ymin = Proportion_lwr, 
                      ymax = Proportion_upr)) +
  facet_grid(~ Treatment, 
             scales = "free", 
             space = "free") +
  theme_minimal()

deco_summary %>%
  ggplot() +
  geom_pointrange(aes(Day, Proportion_mean, 
                      ymin = Proportion_lwr, 
                      ymax = Proportion_upr)) +
  facet_grid(vars(Tank), 
             scales = "free", 
             space = "free") +
  theme_minimal()
# Each tank has enough information to estimate all
# parameters.

# 2.2 Prior simulation ####
# Here I use a model I developed specifically for macroalgal 
# decomposition (github.com/lukaseamus/limbodeco).

# R doesn't have a built-in log1p_exp function
log1p_exp <- function(x) {
  ifelse(
    x > 0, 
    x + log1p(exp(-x)),
    log1p(exp(x))
  )
}

# I am taking the mean k value for Ecklonia radiata from 
# Simpkins et al. 2025 (doi: 10.1002/lno.70006) as my prior
# for tau: 0.06 d^-1. Wright et al. 2024 (doi: 10.1093/aob/mcad167)
# show that the half-life of detrital photosynthesis (mu) ranges from
# 15 to 83 d across treatments, 40 d seems like a reasonable mean. 
# The  relative growth rate of E. radiata is typically expected to be less 
# than 0.01 d^-1 (Fairhead & Cheshire 2004, doi: 10.1007/s00227-004-1308-8) 
# and is often negative under suboptimal condistions (Xiao et al. 2015, doi:
# 10.1371/journal.pone.0143031). But since alpha parameterises the initial
# relative growth rate and all detritus was known to be alive at the
# start of the experiment (Wright et al. 2024, doi: 10.1093/aob/mcad167),
# it makes sense to constrain it to positive values close between
# zero and 0.01. I will judge this based on prior simulation.

# There is no information to inform epsilon because there are no
# data at t = 0 but I won't treat it as constant. Instead I will
# give it a prior and use partial pooling. The initial precision
initial_posterior %>%
  mutate(Initial_p = Initial / mean(Initial)) %>%
  summarise(Initial_p_mean = mean(Initial_p),
            Initial_p_sd = sd(Initial_p),
            Initial_p_nu = Initial_p_mean * ( 1 + Initial_p_mean ) / Initial_p_sd^2,
            n = n())
# should not be used because initial variation is already accounted
# for in the standard deviation of the proportion. Since variability
# divided by variability is zero
initial_posterior %>%
  mutate(Proportion = Initial / Initial) %>%
  summarise(Proportion_mean = mean(Proportion),
            Proportion_sd = sd(Proportion),
            n = n())
# the actual precision should be high. I will go with the default 
# values for epsilon, lambda and tau proposed previously 
# (github.com/lukaseamus/limbodeco).

require(extraDistr)
tibble(n = 1:1e3,
       log_alpha_mu = rnorm( 1e3 , log(0.001) , 0.3 ), 
       log_mu_mu = rnorm( 1e3 , log(40) , 0.3 ),
       log_tau_mu = rnorm( 1e3 , log(0.06) , 0.3 ),
       log_alpha_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ), 
       log_mu_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_tau_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_epsilon_mu = rnorm( 1e3 , log(4e4) , 0.3 ),
       log_lambda_mu = rnorm( 1e3 , log(0.1) , 0.3 ),
       log_theta_mu = rnorm( 1e3 , log(500) , 0.3 ),
       log_epsilon_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_lambda_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_theta_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       # There are two factors in the multilevel model
       alpha = exp( rnorm( 1e3 , log_alpha_mu , log_alpha_sigma ) +
                      rnorm( 1e3 , 0 , log_alpha_sigma ) ),
       mu = exp( rnorm( 1e3 , log_mu_mu , log_mu_sigma ) +
                   rnorm( 1e3 , 0 , log_mu_sigma ) ),
       tau = exp( rnorm( 1e3 , log_tau_mu , log_tau_sigma ) +
                    rnorm( 1e3 , 0 , log_tau_sigma ) ),
       epsilon = exp( rnorm( 1e3 , log_epsilon_mu , log_epsilon_sigma ) +
                        rnorm( 1e3 , 0 , log_epsilon_sigma ) ),
       lambda = exp( rnorm( 1e3 , log_lambda_mu , log_lambda_sigma ) +
                       rnorm( 1e3 , 0 , log_lambda_sigma ) ),
       theta = exp( rnorm( 1e3 , log_theta_mu , log_theta_sigma ) + 
                      rnorm( 1e3 , 0 , log_theta_sigma ) )) %>%
  expand_grid(Day = deco_summary %$% 
                seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    p_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    nu = theta + (epsilon - theta) * exp( -lambda * Day ),
    p = rbetapr( n() , p_mu * ( 1 + nu ) , 2 + nu )
  ) %>%
  pivot_longer(cols = c(p_mu, k, nu, p),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 2.3 Stan model ####
deco_c_model <- here("Decomposition", "Stan", "deco_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

deco_nc_model <- here("Decomposition", "Stan", "deco_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

deco_c_samples <- deco_c_model$sample(
          data = deco_summary %>%
            select(Day, Proportion_mean, Proportion_sd,
                   Treatment, Tank) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

deco_nc_samples <- deco_nc_model$sample(
          data = deco_summary %>%
            select(Day, Proportion_mean, Proportion_sd,
                   Treatment, Tank) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# 2.4 Model checks ####
# Rhat
deco_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 100% of rhat above 1.001. rhat = 1.27 ± 0.0711.

deco_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# Less than 1% of rhat above 1.001. rhat = 1.00 ± 0.000143.
# Already there's a clear winner.

# Chains
deco_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay()

deco_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay()
# Chains are also far better for non-centred model.

# Pairs
deco_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_alpha_t[1]", "log_mu_t[1]", "log_tau_t[1]"))
deco_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_alpha_t[2]", "log_mu_t[2]", "log_tau_t[2]"))
deco_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_alpha_t[3]", "log_mu_t[3]", "log_tau_t[3]"))
deco_c_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_alpha_t[4]", "log_mu_t[4]", "log_tau_t[4]"))
# Posteriors look spiky. Only some weak correlation.

deco_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_alpha_t[1]", "log_mu_t[1]", "log_tau_t[1]"))
deco_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_alpha_t[2]", "log_mu_t[2]", "log_tau_t[2]"))
deco_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_alpha_t[3]", "log_mu_t[3]", "log_tau_t[3]"))
deco_nc_samples$draws(format = "df") %>%
  mcmc_pairs(pars = c("log_alpha_t[4]", "log_mu_t[4]", "log_tau_t[4]"))
# Posteriors look smooth. Barely any correlation.

# 2.5 Prior-posterior comparison ####
# Hierarchical priors cannot effectively sampled hen centred.
# Hence the non-centred model will be used to sample priors.
deco_prior <- prior_samples(
  model = deco_nc_model,
  data = deco_summary %>%
    select(Day, Proportion_mean, Proportion_sd,
           Treatment, Tank) %>%
    compose_data()
)

# Centred model
deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_c_samples,
    group = deco_summary %>% select(Treatment),
    parameters = c("log_alpha_mu", "log_alpha_sigma_t",
                   "log_alpha_t[Treatment]",
                   "log_mu_mu", "log_mu_sigma_t",
                   "log_mu_t[Treatment]",
                   "log_tau_mu", "log_tau_sigma_t",
                   "log_tau_t[Treatment]",
                   "log_epsilon_mu", "log_epsilon_sigma_t",
                   "log_epsilon_t[Treatment]",
                   "log_lambda_mu", "log_lambda_sigma_t",
                   "log_lambda_t[Treatment]",
                   "log_theta_mu", "log_theta_sigma_t",
                   "log_theta_t[Treatment]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Treatment")

deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_c_samples,
    group = deco_summary %>% select(Tank),
    parameters = c("log_alpha_sigma_ta", "log_alpha_ta[Tank]",
                   "log_mu_sigma_ta", "log_mu_ta[Tank]", 
                   "log_tau_sigma_ta", "log_tau_ta[Tank]",
                   "log_epsilon_sigma_ta", "log_epsilon_ta[Tank]",
                   "log_lambda_sigma_ta", "log_lambda_ta[Tank]",
                   "log_theta_sigma_ta", "log_theta_ta[Tank]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Tank", ridges = T)
# Again spiky posteriors.

# Non-centred model
deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_nc_samples,
    group = deco_summary %>% select(Treatment),
    parameters = c("log_alpha_mu", "log_alpha_sigma_t",
                   "log_alpha_t[Treatment]",
                   "log_mu_mu", "log_mu_sigma_t",
                   "log_mu_t[Treatment]",
                   "log_tau_mu", "log_tau_sigma_t",
                   "log_tau_t[Treatment]",
                   "log_epsilon_mu", "log_epsilon_sigma_t",
                   "log_epsilon_t[Treatment]",
                   "log_lambda_mu", "log_lambda_sigma_t",
                   "log_lambda_t[Treatment]",
                   "log_theta_mu", "log_theta_sigma_t",
                   "log_theta_t[Treatment]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Treatment")

deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_nc_samples,
    group = deco_summary %>% select(Tank),
    parameters = c("log_alpha_sigma_ta", "log_alpha_ta[Tank]",
                   "log_mu_sigma_ta", "log_mu_ta[Tank]", 
                   "log_tau_sigma_ta", "log_tau_ta[Tank]",
                   "log_epsilon_sigma_ta", "log_epsilon_ta[Tank]",
                   "log_lambda_sigma_ta", "log_lambda_ta[Tank]",
                   "log_theta_sigma_ta", "log_theta_ta[Tank]"),
    format = "long"
  ) %>%
  prior_posterior_plot(group_name = "Tank", ridges = T)
# Posteriors for the non-centred model are smoother.
# Proceed with that.

# 2.6 Parameters ####
# Global parameters
deco_prior_posterior_global <- deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_nc_samples,
    parameters = c("log_alpha_mu", "log_alpha_sigma_t", "log_alpha_sigma_ta",
                   "log_mu_mu", "log_mu_sigma_t", "log_mu_sigma_ta",
                   "log_tau_mu", "log_tau_sigma_t", "log_tau_sigma_ta",
                   "log_epsilon_mu", "log_epsilon_sigma_t", "log_epsilon_sigma_ta",
                   "log_lambda_mu", "log_lambda_sigma_t", "log_lambda_sigma_ta",
                   "log_theta_mu", "log_theta_sigma_t", "log_theta_sigma_ta"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for unobserved treatments and tanks
    alpha = exp(
      rnorm( n() , log_alpha_mu , log_alpha_sigma_t ) + 
        rnorm( n() , 0 , log_alpha_sigma_ta )
    ),
    mu = exp(
      rnorm( n() , log_mu_mu , log_mu_sigma_t ) +
        rnorm( n() , 0 , log_mu_sigma_ta )
    ),
    tau = exp(
      rnorm( n() , log_tau_mu , log_tau_sigma_t ) +
        rnorm( n() , 0 , log_tau_sigma_ta )
    ),
    epsilon = exp(
      rnorm( n() , log_epsilon_mu , log_epsilon_sigma_t ) +
        rnorm( n() , 0 , log_epsilon_sigma_ta )
    ),
    lambda = exp(
      rnorm( n() , log_lambda_mu , log_lambda_sigma_t ) +
        rnorm( n() , 0 , log_lambda_sigma_ta )
    ),
    theta = exp(
      rnorm( n() , log_theta_mu , log_theta_sigma_t ) +
        rnorm( n() , 0 , log_theta_sigma_ta )
    )
  ) %>%
  select(starts_with("."), distribution, 
         alpha, mu, tau, epsilon, lambda, theta) %T>%
  print()

deco_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  group_by(distribution, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n())

# Treatment parameters
deco_prior_posterior <- deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_nc_samples,
    group = deco_summary %>% select(Treatment),
    parameters = c("log_alpha_sigma_ta", "log_alpha_t[Treatment]",
                   "log_mu_sigma_ta", "log_mu_t[Treatment]",
                   "log_tau_sigma_ta", "log_tau_t[Treatment]",
                   "log_epsilon_sigma_ta", "log_epsilon_t[Treatment]",
                   "log_lambda_sigma_ta", "log_lambda_t[Treatment]",
                   "log_theta_sigma_ta", "log_theta_t[Treatment]"),
    format = "short"
  ) %>% 
  mutate( # Calculate parameters for new tanks
    alpha = exp( log_alpha_t + rnorm( n() , 0 , log_alpha_sigma_ta ) ),
    mu = exp( log_mu_t + rnorm( n() , 0 , log_mu_sigma_ta ) ),
    tau = exp( log_tau_t + rnorm( n() , 0 , log_tau_sigma_ta ) ),
    epsilon = exp( log_epsilon_t + rnorm( n() , 0 , log_epsilon_sigma_ta ) ),
    lambda = exp( log_lambda_t + rnorm( n() , 0 , log_lambda_sigma_ta ) ),
    theta = exp( log_theta_t + rnorm( n() , 0 , log_theta_sigma_ta ) )
  ) %>% # Remove redundant priors
  filter(!(Treatment %>% str_detect("Light") & 
             distribution == "prior")) %>%
  mutate( # Embed prior in treatment
    Treatment = if_else(
      distribution == "prior", "Prior", Treatment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Treatment, 
         alpha, mu, tau, epsilon, lambda, theta) %T>%
  print()

deco_prior_posterior %>%
  pivot_longer(cols = -c(starts_with("."), Treatment),
               names_to = "parameter") %>%
  group_by(Treatment, parameter) %>%
  summarise(mean = mean(value), sd = sd(value), n = n()) %>%
  print(n = 36)

# 2.7 Prediction ####
# Predict across predictor range
deco_prediction <- deco_prior_posterior %>%
  spread_continuous(data = deco_summary %>%
                      # Ensure predictor range starts at 0
                      bind_rows(
                        expand_grid(
                          Day = 0, 
                          Proportion_mean = 1,
                          Treatment = deco_summary %$% 
                            levels(Treatment) %>% fct()
                        )
                      ),
                    group_name = "Treatment",
                    predictor_name = "Day", 
                    length = 200) %>%
  mutate(
    p_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    # Half-life only makes sense for negative k values (decomposition).
    t0.5 = if_else( k < 0 , log(0.5) / k , Inf ),  
    nu = ( epsilon - theta ) * exp( -lambda * Day ) + theta,
    p = rbetapr( n() , p_mu * ( 1 + nu ) , 2 + nu )
  ) %T>%
  print()

# Summarise predictions
deco_prediction_summary <- deco_prediction %>%
  group_by(Day, Treatment) %>%
  median_qi(p_mu, k, nu, p, .width = c(.5, .8, .9)) %T>%
  print()

# Half-life needs to be summarised separately because of infinities
t0.5_summary <- deco_prediction %>%
  filter(is.finite(t0.5)) %>% # drop Inf
  group_by(Day, Treatment) %>%
  median_qi(t0.5, .width = c(.5, .8, .9)) %T>%
  print()

deco_prediction_summary %<>%
  left_join(
    t0.5_summary,
    by = c("Day", "Treatment", ".width", 
           ".point", ".interval")
  ) %T>%
  print()

# 3. Visualise data ####
# Define custom theme
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = margin(0.2, 0.5, 0.2, 0.2, unit = "cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 12, hjust = 0),
                 axis.text = element_text(size = 10, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black", lineend = "square"),
                 legend.key = element_blank(),
                 legend.key.width = unit(.25, "cm"),
                 legend.key.height = unit(.45, "cm"),
                 legend.key.spacing.x = unit(.5, "cm"),
                 legend.key.spacing.y = unit(.05, "cm"),
                 legend.background = element_blank(),
                 legend.position = "top",
                 legend.justification = 0,
                 legend.text = element_text(size = 12, hjust = 0),
                 legend.title = element_blank(),
                 legend.margin = margin(0, 0, 0, 0, unit = "cm"),
                 strip.background = element_blank(),
                 strip.text = element_text(size = 12, hjust = 0),
                 panel.spacing = unit(1, "cm"),
                 text = element_text(family = "Futura"))

# Annotation
require(glue)
deco_annotation <- deco_summary %>%
  group_by(Treatment) %>%
  summarise(n = glue("italic(n)*' = {n()}'")) %T>%
  print()

# Plot
Fig_1a <- deco_prediction_summary %>%
  filter(Treatment != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = 100) +
    geom_pointrange(data = deco_summary,
                    aes(Day, Proportion_mean * 100,
                        ymin = Proportion_lwr * 100,
                        ymax = Proportion_upr * 100,
                        colour = Treatment),
                    shape = 16, alpha = 0.5) +
    geom_line(aes(Day, p * 100, colour = Treatment)) +
    geom_ribbon(aes(Day, ymin = p.lower * 100, ymax = p.upper * 100,
                    alpha = factor(.width), fill = Treatment)) +
    geom_text(data = deco_annotation,
              aes(c(120, 120, 60, 60), 170, label = n),
              family = "Futura", parse = TRUE,
              size.unit = "pt", size = 10, hjust = 1) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = seq(0, 120, 30)) +
    scale_y_continuous(breaks = seq(0, 180, 60)) +
    labs(x = "Detrital age (days)",
         y = "Dry mass (%)") +
    coord_cartesian(ylim = c(0, 180), expand = F, clip = "off") +
    facet_grid(~ Treatment, space = "free_x", scales = "free_x") +
    mytheme + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())

Fig_1b <- deco_prediction_summary %>%
  filter(Treatment != "Prior") %>%
  ggplot() +
    geom_line(aes(Day, k, colour = Treatment)) +
    geom_ribbon(aes(Day, ymin = k.lower, ymax = k.upper,
                    alpha = factor(.width), fill = Treatment)) +
    # I include this text out of bounds to maintain x scale limits.
    geom_text(data = deco_annotation, 
              aes(c(120, 120, 60, 60), 1, label = n),
              family = "Futura", parse = TRUE,
              size.unit = "pt", size = 10, hjust = 1) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = seq(0, 120, 30)) +
    scale_y_continuous(breaks = seq(-0.12, 0, 0.04),
                       labels = scales::label_number(
                         accuracy = c(rep(0.01, 3), 1),
                         style_negative = "minus"
                         )
                       ) +
    labs(x = "Detrital age (days)",
         y = expression(italic(k)*" (d"^-1*")")) +
    coord_cartesian(ylim = c(-0.12, 0), expand = F, clip = "off") +
    facet_grid(~ Treatment, space = "free_x", scales = "free_x") +
    mytheme + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          strip.text = element_blank())

Fig_1c <- deco_prediction_summary %>%
  filter(Treatment != "Prior") %>%
  ggplot() +
    geom_line(aes(Day, t0.5, colour = Treatment)) +
    geom_ribbon(aes(Day, ymin = .lower, ymax = .upper,
                    alpha = factor(.width), fill = Treatment)) +
    # I include this text out of bounds to maintain x scale limits.
    geom_text(data = deco_annotation, 
              aes(c(120, 120, 60, 60), 1000, label = n),
              family = "Futura", parse = TRUE,
              size.unit = "pt", size = 10, hjust = 1) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_x_continuous(breaks = seq(0, 120, 30)) +
    scale_y_continuous(breaks = seq(0, 120, 40)) +
    labs(x = "Detrital age (days)",
         y = "Half-life (days)") +
    coord_cartesian(ylim = c(0, 120), xlim = c(0, NA), expand = F) +
    facet_grid(~ Treatment, space = "free_x", scales = "free_x") +
    mytheme +
    theme(strip.text = element_blank())

require(patchwork)
Fig_1 <- ( Fig_1a  / Fig_1b  / Fig_1c ) +
  plot_annotation(tag_levels = c("a", "b", "c")) &
  theme(plot.tag = element_text(family = "Futura", size = 15, face = "bold"))

Fig_1 %>%
  ggsave(filename = "Fig_1.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 20, units = "cm")
