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

# 1.2 Estimate initial dry mass ####
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

# R-hat
initial_samples$summary() %>%
  summarise(rhat_1.001 = mean( rhat > 1.001 ),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))

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

# Compare to data summary
deco %>%
  filter(Day == 0) %>%
  summarise(Initial_mean = mean(Dry),
            Initial_sd = sd(Dry),
            n = n())
# Mean is similar but variance is underestimated

# 1.3 Calculations ####
deco %<>%
  filter(Day != 0) %>%
  cross_join(
    initial_posterior %>% 
      select(starts_with("."), Initial)
  ) %T>%
  print()

deco %<>%
  mutate(Ratio = Dry / Initial,
         k = -log(Ratio) / Day,
         dm_dt = -k * Dry)

deco_summary <- deco %>%
  summarise(
    across(
      c(Initial, Ratio, k, dm_dt),
      list(
        mean = mean, sd = sd,
        lower = ~qi(.x, .width = 0.90)[1],
        upper = ~qi(.x, .width = 0.90)[2])
    ),
    n = n(),
    .by = c(Date, Tank, Temperature, PAR, Wet, Dry, d15N, d13C,
            N, C, C_N, Day, Treatment)
  ) %T>%
  print(n = 122)

deco_summary %>% # Save
  write_rds(here("Decomposition", "RDS", "deco_summary.rds"))

# 2. Macroalgal model ####
# 2.1 Data exploration ####
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

deco_summary %>%
  ggplot() +
  geom_pointrange(aes(Day, Ratio_mean, 
                      ymin = Ratio_lower, 
                      ymax = Ratio_upper)) +
  facet_grid(~ Treatment, 
             scales = "free", 
             space = "free") +
  mytheme

deco_summary %>%
  ggplot() +
  geom_pointrange(aes(Day, Ratio_mean, 
                      ymin = Ratio_lower, 
                      ymax = Ratio_upper)) +
  facet_wrap(~Tank) +
  mytheme
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
# 15 to 83 d across treatments, 50 d seems like a reasonable mean. 
# The  relative growth rate of E. radiata is typically expected to be less 
# than 0.01 d^-1 (Fairhead & Cheshire 2004, doi: 10.1007/s00227-004-1308-8) 
# and is often negative under suboptimal condistions (Xiao et al. 2015, doi:
# 10.1371/journal.pone.0143031). I will parameterise in terms of
# delta = alpha + tau, and set delta < tau, so prior median alpha is negative.

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
  mutate(Ratio = Initial / Initial) %>%
  summarise(Ratio_mean = mean(Ratio),
            Ratio_sd = sd(Ratio),
            n = n())
# the actual precision should be high. I will go with the default 
# values for epsilon, lambda and tau proposed previously 
# (github.com/lukaseamus/limbodeco).

require(extraDistr)
tibble(n = 1:1e3,
       log_delta_mu = rnorm( 1e3 , log(0.05) , 0.3 ), # delta = alpha + tau
       log_mu_mu = rnorm( 1e3 , log(50) , 0.3 ),
       log_tau_mu = rnorm( 1e3 , log(0.06) , 0.3 ),
       log_delta_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ), 
       log_mu_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       log_tau_sigma = rtnorm( 1e3 , 0 , 0.3 , 0 ),
       # There are two factors in the multilevel model (treatment and tank)
       delta = exp( rnorm( 1e3 , log_delta_mu , log_delta_sigma ) +
                      rnorm( 1e3 , 0 , log_delta_sigma ) ),
       mu = exp( rnorm( 1e3 , log_mu_mu , log_mu_sigma ) +
                   rnorm( 1e3 , 0 , log_mu_sigma ) ),
       tau = exp( rnorm( 1e3 , log_tau_mu , log_tau_sigma ) +
                    rnorm( 1e3 , 0 , log_tau_sigma ) ),
       alpha = delta - tau,
       epsilon = rgamma( 1e3 , 4e4^2 / 2e4^2 , 4e4 / 2e4^2 ),
       lambda = rexp( 1e3 , 1 ),
       theta = rgamma( 1e3 , 500^2 / 250^2 , 500 / 250^2 )) %>%
  expand_grid(Day = deco_summary %$% 
                seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    nu = theta + (epsilon - theta) * exp( -lambda * Day ),
    r = rbetapr( n() , r_mu * ( 1 + nu ) , 2 + nu )
  ) %>%
  pivot_longer(cols = c(r_mu, k, nu, r),
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
            select(Day, Ratio_mean, Ratio_sd,
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
            select(Day, Ratio_mean, Ratio_sd,
                   Treatment, Tank) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# Save draws
deco_c_samples$draws() %>%
  write_rds(here("Decomposition", "RDS", "deco_c_samples.rds"))
deco_c_samples$draws(format = "df") %>%
  write_rds(here("Decomposition", "RDS", "deco_c_samples_df.rds"))

deco_nc_samples$draws() %>%
  write_rds(here("Decomposition", "RDS", "deco_nc_samples.rds"))
deco_nc_samples$draws(format = "df") %>%
  write_rds(here("Decomposition", "RDS", "deco_nc_samples_df.rds"))

# 2.4 Model checks ####
# Rhat
deco_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 100% of rhat above 1.001. rhat = 1.25 ± 0.0622.

deco_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.000155.
# Already there's a clear winner.

# Chains
deco_c_chains <- deco_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off") +
  mytheme

deco_c_chains %>%
  ggsave(filename = "deco_c_chains.pdf", path = here("Decomposition", "Plots"),
         device = cairo_pdf, width = 100*14/16, height = 100*14/16, units = "cm")

deco_nc_chains <- deco_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Non-centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off") +
  mytheme

deco_nc_chains %>%
  ggsave(filename = "deco_nc_chains.pdf", path = here("Decomposition", "Plots"),
         device = cairo_pdf, width = 100, height = 100, units = "cm")
# Chains are much better for non-centred model

rm(deco_c_chains, deco_nc_chains) # Clean up
gc()

# Pairs
deco_c_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("log_delta_mu", "log_delta_sigma_t", "log_delta_t[1]", "log_delta_t[2]",
             "log_delta_sigma_ta", "log_delta_ta[2]", "log_delta_ta[10]",
             "log_mu_mu", "log_mu_sigma_t", "log_mu_t[1]", "log_mu_t[2]",
             "log_mu_sigma_ta", "log_mu_ta[2]", "log_mu_ta[10]",
             "log_tau_mu", "log_tau_sigma_t", "log_tau_t[1]", "log_tau_t[2]",
             "log_tau_sigma_ta", "log_tau_ta[2]", "log_tau_ta[10]",
             "epsilon", "lambda", "theta"),
    grid_args = list(top = "Centred model")
  ) %>%
  ggsave(filename = "deco_c_pairs.png", path = here("Decomposition", "Plots"),
         width = 100, height = 100, units = "cm", bg = "white")

deco_nc_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("log_delta_mu", "log_delta_sigma_t", "log_delta_t[1]", "log_delta_t[2]",
             "log_delta_sigma_ta", "log_delta_ta[2]", "log_delta_ta[10]",
             "log_mu_mu", "log_mu_sigma_t", "log_mu_t[1]", "log_mu_t[2]",
             "log_mu_sigma_ta", "log_mu_ta[2]", "log_mu_ta[10]",
             "log_tau_mu", "log_tau_sigma_t", "log_tau_t[1]", "log_tau_t[2]",
             "log_tau_sigma_ta", "log_tau_ta[2]", "log_tau_ta[10]",
             "epsilon", "lambda", "theta"),
    grid_args = list(top = "Non-centred model")
  ) %>%
  ggsave(filename = "deco_nc_pairs.png", path = here("Decomposition", "Plots"),
         width = 100, height = 100, units = "cm", bg = "white")

# 2.5 Prior-posterior comparison ####
# Hierarchical priors cannot effectively sampled hen centred.
# Hence the non-centred model will be used to sample priors.
deco_prior <- prior_samples(
  model = deco_nc_model,
  data = deco_summary %>%
    select(Day, Ratio_mean, Ratio_sd,
           Treatment, Tank) %>%
    compose_data()
)

# Centred model
deco_c_prior_posterior_treatment <- deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_c_samples,
    group = deco_summary %>% select(Treatment, Tank),
    parameters = c("log_delta_mu", "log_delta_sigma_t", "log_delta_t[Treatment]",
                   "log_delta_sigma_ta", "log_mu_mu", "log_mu_sigma_t", "log_mu_t[Treatment]",
                   "log_mu_sigma_ta", "log_tau_mu", "log_tau_sigma_t", "log_tau_t[Treatment]",
                   "log_tau_sigma_ta", "epsilon", "lambda", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Treatment") +
  scale_x_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(title = "Centred model") +
  coord_cartesian(expand = FALSE) +
  mytheme +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

deco_c_prior_posterior_tank <- deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_c_samples,
    group = deco_summary %>% select(Tank),
    parameters = c("log_delta_ta[Tank]", "log_mu_ta[Tank]", "log_tau_ta[Tank]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Tank", ridges = TRUE) +
  scale_x_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  coord_cartesian(expand = FALSE) +
  mytheme +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

# Non-centred model
deco_nc_prior_posterior_treatment <- deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_nc_samples,
    group = deco_summary %>% select(Treatment, Tank),
    parameters = c("log_delta_mu", "log_delta_sigma_t", "log_delta_t[Treatment]",
                   "log_delta_sigma_ta", "log_mu_mu", "log_mu_sigma_t", "log_mu_t[Treatment]",
                   "log_mu_sigma_ta", "log_tau_mu", "log_tau_sigma_t", "log_tau_t[Treatment]",
                   "log_tau_sigma_ta", "epsilon", "lambda", "theta"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Treatment") +
  scale_x_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(title = "Non-centred model") +
  coord_cartesian(expand = FALSE) +
  mytheme +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

deco_nc_prior_posterior_tank <- deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_nc_samples,
    group = deco_summary %>% select(Tank),
    parameters = c("log_delta_ta[Tank]", "log_mu_ta[Tank]", "log_tau_ta[Tank]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Tank", ridges = TRUE) +
  scale_x_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  coord_cartesian(expand = FALSE) +
  mytheme +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

# Combine
require(patchwork)
deco_prior_posterior_plot <- 
  ( deco_c_prior_posterior_treatment / deco_c_prior_posterior_tank +
        plot_layout(heights = c(1, 0.5)) ) | 
  ( deco_nc_prior_posterior_treatment / deco_nc_prior_posterior_tank ) +
        plot_layout(heights = c(1, 0.5))

deco_prior_posterior_plot %>%
  ggsave(filename = "deco_prior_posterior.pdf", 
         path = here("Decomposition", "Plots"),
         device = cairo_pdf, width = 60, height = 40, units = "cm")
# Non-centred model is smoother. Choose as optimal model.

rm(deco_nc_prior_posterior_treatment, deco_nc_prior_posterior_tank,
   deco_c_prior_posterior_treatment, deco_c_prior_posterior_tank,
   deco_prior_posterior_plot, deco_c_model, deco_c_samples) # Clean up
gc()

# 2.6 Parameters ####
# Global parameters
deco_prior_posterior_global <- deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_nc_samples,
    parameters = c("log_delta_mu", "log_delta_sigma_t", "log_delta_sigma_ta",
                   "log_mu_mu", "log_mu_sigma_t", "log_mu_sigma_ta",
                   "log_tau_mu", "log_tau_sigma_t", "log_tau_sigma_ta",
                   "epsilon", "lambda", "theta"),
    format = "short"
  ) %>%
  mutate( # Calculate parameters for unobserved treatments and tanks
    delta = exp(
      rnorm( n() , log_delta_mu , log_delta_sigma_t ) + 
        rnorm( n() , 0 , log_delta_sigma_ta )
    ),
    mu = exp(
      rnorm( n() , log_mu_mu , log_mu_sigma_t ) +
        rnorm( n() , 0 , log_mu_sigma_ta )
    ),
    tau = exp(
      rnorm( n() , log_tau_mu , log_tau_sigma_t ) +
        rnorm( n() , 0 , log_tau_sigma_ta )
    ),
    alpha = delta - tau
  ) %>%
  select(starts_with("."), distribution, 
         alpha, mu, tau, epsilon, lambda, theta) %T>%
  print()

deco_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  summarise(mean = mean(value), sd = sd(value), n = n(),
            .by = c(distribution, parameter))

# Treatment parameters
deco_prior_posterior <- deco_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_nc_samples,
    group = deco_summary %>% select(Treatment),
    parameters = c("log_delta_sigma_ta", "log_delta_t[Treatment]",
                   "log_mu_sigma_ta", "log_mu_t[Treatment]",
                   "log_tau_sigma_ta", "log_tau_t[Treatment]",
                   "epsilon", "lambda", "theta"),
    format = "short"
  ) %>% 
  mutate( # Calculate parameters for new tanks
    delta = exp( rnorm( n() , log_delta_t , log_delta_sigma_ta ) ),
    mu = exp( rnorm( n() , log_mu_t , log_mu_sigma_ta ) ),
    tau = exp( rnorm( n() , log_tau_t , log_tau_sigma_ta ) ),
    alpha = delta - tau
  ) %>% # Remove redundant priors (keep only dark prior)
  filter(!(Treatment %>% str_detect("Light") & distribution == "prior")) %>%
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
  summarise(mean = mean(value), sd = sd(value), n = n(),
            .by = c(Treatment, parameter)) %>%
  print(n = 30)

# Save parameter distributions
deco_prior_posterior_global %>%
  write_rds(here("Decomposition", "RDS", "deco_prior_posterior_global.rds"))
deco_prior_posterior %>%
  write_rds(here("Decomposition", "RDS", "deco_prior_posterior.rds"))

# 2.7 Prediction ####
# Predict across predictor range
deco_prediction <- deco_prior_posterior %>%
  spread_continuous(data = deco_summary %>%
                      add_row(
                        Day = 0, # Ensure predictor range starts at 0
                        Treatment = deco_summary %$% 
                          levels(Treatment) %>% fct()
                      ), 
                    group_name = "Treatment",
                    predictor_name = "Day") %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    k = ( alpha + tau ) / ( 1 + exp( 5 / mu * ( Day - mu ) ) ) - tau,
    dm_dt = k * r_mu,
    # Half-life only makes sense for negative k values (decomposition).
    t0.5 = if_else( k < 0 , log(0.5) / k , Inf ),  
    nu = ( epsilon - theta ) * exp( -lambda * Day ) + theta,
    r = rbetapr( n() , r_mu * ( 1 + nu ) , 2 + nu )
  ) %T>%
  print()

# Summarise
deco_prediction_summary <- deco_prediction %>%
  summarise(
    across(
      c(r_mu, k, dm_dt, nu, r),
      list(median = median, 
           lower_0.5 = ~ qi(.x, .width = .5)[1],
           upper_0.5 = ~ qi(.x, .width = .5)[2],
           lower_0.8 = ~ qi(.x, .width = .8)[1],
           upper_0.8 = ~ qi(.x, .width = .8)[2],
           lower_0.9 = ~ qi(.x, .width = .9)[1],
           upper_0.9 = ~ qi(.x, .width = .9)[2]),
      .names = "{.col}.{.fn}"
    ),
    .by = c(Day, Treatment)
  ) %>%
  rename(r_mu = r_mu.median, k = k.median, dm_dt = dm_dt.median,
         nu = nu.median, r = r.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

# Half-life needs to be summarised separately because of infinities
t0.5_prediction_summary <- deco_prediction %>%
  filter(is.finite(t0.5)) %>% # drop Inf
  summarise(
    across(
      t0.5,
      list(median = median, 
           lower_0.5 = ~ qi(.x, .width = .5)[1],
           upper_0.5 = ~ qi(.x, .width = .5)[2],
           lower_0.8 = ~ qi(.x, .width = .8)[1],
           upper_0.8 = ~ qi(.x, .width = .8)[2],
           lower_0.9 = ~ qi(.x, .width = .9)[1],
           upper_0.9 = ~ qi(.x, .width = .9)[2]),
      .names = "{.col}.{.fn}"
    ),
    .by = c(Day, Treatment)
  ) %>%
  rename(t0.5 = t0.5.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

# Clean up
rm(deco_prediction)

# Join summaries
deco_prediction <- deco_prediction_summary %>%
  full_join(t0.5_prediction_summary) %T>%
  print()

# Save predictions
deco_prediction %>%
  write_rds(here("Decomposition", "RDS", "deco_prediction.rds"))



# 3. Conventional model ####
# 3.1 Prior simulation ####
# Here I use the conventional model exp(-k * Day).
# I am taking the mean k value for Ecklonia radiata from 
# Simpkins et al. 2025 (doi: 10.1002/lno.70006) as my prior
# for tau: 0.06 d^-1.

tibble(n = 1:1e3,
       log_k_mu = rnorm( 1e3 , log(0.06) , 1 ),
       log_k_sigma = rtnorm( 1e3 , 0 , 1 , 0 ),
       k = exp( rnorm( 1e3 , log_k_mu , log_k_sigma ) +
                  rnorm( 1e3 , 0 , log_k_sigma ) ),
       sigma = rexp( 1e3 , 1 )) %>%
  expand_grid(Day = deco_summary %$% 
                seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    r_mu = exp(-k * Day),
    r = rnorm( n() , r_mu , sigma )
  ) %>%
  pivot_longer(cols = c(r_mu, r),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 3.2 Stan model ####
deco_k_model <- here("Decomposition", "Stan", "deco_k.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

deco_k_samples <- deco_k_model$sample(
          data = deco_summary %>%
            select(Day, Ratio_mean, Treatment, Tank) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print(max_rows = 200)

# Save draws
deco_k_samples$draws() %>%
  write_rds(here("Decomposition", "RDS", "deco_k_samples.rds"))
deco_k_samples$draws(format = "df") %>%
  write_rds(here("Decomposition", "RDS", "deco_k_samples_df.rds"))

# 3.3 Model checks ####
# Rhat
deco_k_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000718.

# Chains
deco_k_chains <- deco_k_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Conventional model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off") +
  mytheme

deco_k_chains %>%
  ggsave(filename = "deco_k_chains.pdf", path = here("Decomposition", "Plots"),
         device = cairo_pdf, width = 40, height = 40, units = "cm")

rm(deco_k_chains) # Clean up
gc()

# Pairs
deco_k_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("log_k_mu", "log_k_sigma_t", "log_k_t[1]", "log_k_t[2]",
             "log_k_sigma_ta", "log_k_ta[2]", "log_k_ta[10]", "sigma"),
    grid_args = list(top = "Conventional model")
  ) %>%
  ggsave(filename = "deco_k_pairs.png", path = here("Decomposition", "Plots"),
         width = 40, height = 40, units = "cm", bg = "white")

# 3.4 Prior-posterior comparison ####
deco_k_prior <- prior_samples(
  model = deco_k_model,
  data = deco_summary %>%
    select(Day, Ratio_mean, Treatment, Tank) %>%
    compose_data()
)

deco_k_prior_posterior_treatment <- deco_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_k_samples,
    group = deco_summary %>% select(Treatment, Tank),
    parameters = c("log_k_mu", "log_k_sigma_t", "log_k_t[Treatment]",
                   "log_k_sigma_ta", "sigma"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Treatment") +
  scale_x_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  labs(title = "Conventional model") +
  coord_cartesian(expand = FALSE) +
  mytheme +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

deco_k_prior_posterior_tank <- deco_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_k_samples,
    group = deco_summary %>% select(Tank),
    parameters = c("log_k_ta[Tank]"),
    format = "long"
    ) %>%
  prior_posterior_plot(group_name = "Tank", ridges = TRUE) +
  scale_x_continuous(
    labels = scales::label_number(style_negative = "minus")
  ) +
  coord_cartesian(expand = FALSE) +
  mytheme +
  theme(axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_blank())

# Combine
deco_k_prior_posterior_plot <- 
  ( deco_k_prior_posterior_treatment / deco_k_prior_posterior_tank +
        plot_layout(heights = c(1, 0.5)) )

deco_k_prior_posterior_plot %>%
  ggsave(filename = "deco_k_prior_posterior.pdf", 
         path = here("Decomposition", "Plots"),
         device = cairo_pdf, width = 20, height = 25, units = "cm")

rm(deco_k_prior_posterior_treatment, deco_k_prior_posterior_tank,
   deco_k_prior_posterior_plot) # Clean up
gc()

# 3.5 Parameters ####
# Global parameters
deco_k_prior_posterior_global <- deco_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_k_samples,
    parameters = c("log_k_mu", "log_k_sigma_t", "log_k_sigma_ta", "sigma"),
    format = "short"
  ) %>%
  mutate( # Calculate k for unobserved treatments and tanks
    k = exp(
      rnorm( n() , log_k_mu , log_k_sigma_t ) +
        rnorm( n() , 0 , log_k_sigma_ta )
    )
  ) %>%
  select(starts_with("."), distribution, k, sigma) %T>%
  print()

deco_k_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  summarise(mean = mean(value), sd = sd(value), n = n(),
            .by = c(distribution, parameter))

# Treatment parameters
deco_k_prior_posterior <- deco_k_prior %>% 
  prior_posterior_draws(
    posterior_samples = deco_k_samples,
    group = deco_summary %>% select(Treatment),
    parameters = c("log_k_sigma_ta", "log_k_t[Treatment]", "sigma"),
    format = "short"
  ) %>% 
  mutate( # Calculate k for new tanks
    k = exp( rnorm( n() , log_k_t , log_k_sigma_ta ) )
  ) %>% # Remove redundant priors (keep only dark prior)
  filter(!(Treatment %>% str_detect("Light") & distribution == "prior")) %>%
  mutate( # Embed prior in treatment
    Treatment = if_else(
      distribution == "prior", "Prior", Treatment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Treatment, k, sigma) %T>%
  print()

deco_k_prior_posterior %>%
  pivot_longer(cols = -c(starts_with("."), Treatment),
               names_to = "parameter") %>%
  summarise(mean = mean(value), sd = sd(value), n = n(),
            .by = c(Treatment, parameter)) %>%
  print(n = 30)

# Save parameter distributions
deco_k_prior_posterior_global %>%
  write_rds(here("Decomposition", "RDS", "deco_k_prior_posterior_global.rds"))
deco_k_prior_posterior %>%
  write_rds(here("Decomposition", "RDS", "deco_k_prior_posterior.rds"))

# 3. Tables ####
# Join prior_posterior
deco_prior_posterior_joined <- deco_prior_posterior %>%
  full_join(deco_k_prior_posterior) %>%
  # Calculate half-life of dry mass based on k and tau (final k)
  mutate(t0.5_k = log(2)/k,
         t0.5_tau = log(2)/tau,
         k_diff = tau - k,
         t0.5_diff = t0.5_k - t0.5_tau) %T>%
  print()

require(glue)
Table_1 <- deco_prior_posterior_joined %>%
  mutate(
    k = k * 100, # Convert exponential rates to %
    tau = tau * 100
  ) %>%
  summarise(
    across(
      c(alpha, k, tau, t0.5_k, t0.5_tau, k_diff, t0.5_diff),
      list(mean = mean, sd = sd, median = median)
    ),
    P = mean( k_diff > 0 ),
    P_alpha = mean( alpha > 0 ),
    n = n(),
    .by = c(Treatment)
  ) %>%
  mutate(
    across(where(is.numeric), ~signif(.x, 2)),
    k = glue("{k_mean} ± {k_sd} ({k_median})"),
    tau = glue("{tau_mean} ± {tau_sd} ({tau_median})"),
    t0.5_k = glue("{t0.5_k_mean} ± {t0.5_k_sd} ({t0.5_k_median})"),
    t0.5_tau = glue("{t0.5_tau_mean} ± {t0.5_tau_sd} ({t0.5_tau_median})"),
    t0.5_diff = glue("{t0.5_diff_mean} ± {t0.5_diff_sd} ({t0.5_diff_median})")
  ) %>%
  select(!c(ends_with("mean"), ends_with("median"), ends_with("sd"))) %T>%
  print()

Table_1 %>%
  write_csv(here("Tables", "Table_1.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = Table_1) %>%
  print(target = here("Tables", "Table_1.docx"))

deco_contrast <- deco_prior_posterior_joined %>%
  filter(Treatment != "Prior") %>%
  select(starts_with("."), Treatment, alpha, tau, k, t0.5_tau, t0.5_k) %>%
  mutate(alpha = alpha * 100,
         tau = tau * 100,
         k = k * 100) %>%
  pivot_longer(cols = c(alpha, tau, k, t0.5_tau, t0.5_k),
               names_to = "parameter") %>%
  pivot_wider(names_from = Treatment,
              values_from = value) %>%
  # I want differences relative to ideal (light 15°C)
  mutate(D15vL15_diff = `Dark 15°C` - `Light 15°C`,
         D15vL15_ratio = `Dark 15°C` / `Light 15°C`,
         L20vL15_diff = `Light 20°C` - `Light 15°C`,
         L20vL15_ratio = `Light 20°C` / `Light 15°C`,
         L25vL15_diff = `Light 25°C` - `Light 15°C`,
         L25vL15_ratio = `Light 25°C` / `Light 15°C`) %>%
  select(-c(`Dark 15°C`, `Light 15°C`, `Light 20°C`, `Light 25°C`)) %>%
  pivot_longer(cols = c(ends_with("diff"), ends_with("ratio"))) %>%
  # This step takes long:
  separate(name, into = c("contrast", "type"), sep = "_") %>%
  pivot_wider(values_from = value,
              names_from = type) %T>%
  print()

deco_contrast_summary <- deco_contrast %>%
  summarise(
    across(
      c(diff, ratio),
      list(mean = mean, sd = sd, median = median)
    ),
    P = max( mean( diff > 0 ) , mean( diff < 0 ) ),
    n = n(),
    .by = c(parameter, contrast)
  ) %>%
  mutate(
    across(where(is.numeric), ~signif(.x, 2)),
    diff = glue("{diff_mean} ± {diff_sd} ({diff_median})"),
    ratio = glue("{ratio_mean} ± {ratio_sd} ({ratio_median})")
  ) %>%
  select(!c(ends_with("mean"), ends_with("median"), ends_with("sd"))) %T>%
  print()

deco_contrast_summary %>%
  write_csv(here("Tables", "deco_contrast.csv"))

read_docx() %>%
  body_add_table(value = deco_contrast_summary) %>%
  print(target = here("Tables", "deco_contrast.docx"))

# 4. Figures ####
# Annotation

# deco_annotation <- deco_summary %>%
#   group_by(Treatment) %>%
#   summarise(n = glue("italic(n)*' = {n()}'")) %T>%
#   print()

# Plot
require(ggh4x)
Fig_1a <- deco_prediction %>%
  filter(Treatment != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = 100) +
    geom_pointrange(data = deco_summary,
                    aes(Day, Ratio_mean * 100,
                        ymin = Ratio_lower * 100,
                        ymax = Ratio_upper * 100,
                        colour = Treatment),
                    shape = 16, alpha = 0.5) +
    geom_line(aes(Day, r_mu * 100, colour = Treatment)) +
    geom_ribbon(aes(Day, ymin = r_mu.lower * 100, ymax = r_mu.upper * 100,
                    alpha = factor(.width), fill = Treatment)) +
    # geom_ribbon(data = . %>% filter(.width == 0.9),
    #             aes(Day, ymin = r.lower * 100, 
    #                 ymax = r.upper * 100, fill = Treatment), 
    #             alpha = 0.2) +
    # geom_text(data = deco_annotation,
    #           aes(c(120, 120, 60, 60), 170, label = n),
    #           family = "Futura", parse = TRUE,
    #           size.unit = "pt", size = 10, hjust = 1) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                        guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                      guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(breaks = seq(0, 160, 40)) +
    facet_grid(~ Treatment, space = "free", scales = "free") +
    facetted_pos_scales(
      x = list(
        Treatment %in% c("Dark 15°C", "Light 15°C") ~
          scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)),
        Treatment %in% c("Light 20°C", "Light 25°C") ~
          scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 30))
      )
    ) +
    labs(x = "Detrital age (days)", y = "Dry mass (%)") +
    coord_cartesian(ylim = c(0, 160), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_text(margin = margin(r = -0.3, unit = "cm")))

Fig_1a

Fig_1b <- deco_prediction %>%
  filter(Treatment != "Prior") %>%
  ggplot() +
    geom_hline(yintercept = 0) +
    geom_pointrange(data = deco_summary %>% filter(k_mean <= 0.15), # filter out outliers
                    aes(Day, k_mean * 100,
                        ymin = k_lower * 100,
                        ymax = k_upper * 100,
                        colour = Treatment),
                    shape = 16, alpha = 0.5) +
    geom_line(aes(Day, -k * 100, colour = Treatment)) +
    geom_ribbon(aes(Day, ymin = -k.upper * 100, ymax = -k.lower * 100,
                    alpha = factor(.width), fill = Treatment)) +
    geom_line(
      data = deco_k_prior_posterior %>%
        filter(Treatment != "Prior") %>%
        group_by(Treatment) %>%
        median_qi(k, .width = c(.5, .8, .9)) %>%
        mutate(
          Day = if_else(
            Treatment %in% c("Dark 15°C", "Light 15°C"),
            list(c(0, 120)), list(c(0, 60))
          )
        ) %>%
        unnest(Day),
      aes(Day, k * 100, colour = Treatment), linetype = 5
    ) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                        guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                      guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(labels = scales::label_number(style_negative = "minus")) +
    facet_grid(~ Treatment, space = "free", scales = "free") +
    facetted_pos_scales(
      x = list(
        Treatment %in% c("Dark 15°C", "Light 15°C") ~
          scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)),
        Treatment %in% c("Light 20°C", "Light 25°C") ~
          scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 30))
      )
    ) +
    labs(x = "Detrital age (days)", y = expression("Decay (% day"^-1*")")) +
    coord_cartesian(ylim = c(-5, 15), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_text(margin = margin(r = -0.3, unit = "cm"),
                                      vjust = 1.68))

Fig_1b

# Fig_1c <- deco_prediction %>%
#   filter(Treatment != "Prior") %>%
#   ggplot() +
#     geom_line(aes(Day, t0.5, colour = Treatment)) +
#     geom_ribbon(aes(Day, ymin = t0.5.lower, ymax = t0.5.upper,
#                     alpha = factor(.width), fill = Treatment)) +
#     scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
#                         guide = "none") +
#     scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
#                       guide = "none") +
#     scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
#     scale_x_continuous(breaks = seq(0, 120, 30)) +
#     scale_y_log10() +
#     facet_grid(~ Treatment, space = "free", scales = "free") +
#     labs(x = "Detrital age (days)", y = "Half-life (days)") +
#     coord_cartesian(expand = F) +
#     mytheme
# 
# Fig_1c

# Fig_1c <- deco_prediction %>%
#   filter(Treatment != "Prior") %>%
#   ggplot() +
#     geom_hline(yintercept = 0) +
#     geom_pointrange(data = deco_summary,
#                     aes(Day, dm_dt_mean,
#                         ymin = dm_dt_lower,
#                         ymax = dm_dt_upper,
#                         colour = Treatment),
#                     shape = 16, alpha = 0.5) +
#     geom_line(aes(Day, dm_dt, colour = Treatment)) +
#     geom_ribbon(aes(Day, ymin = dm_dt.lower, ymax = dm_dt.upper,
#                     alpha = factor(.width), fill = Treatment)) +
#     scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
#                         guide = "none") +
#     scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
#                       guide = "none") +
#     scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
#     scale_y_continuous(labels = scales::label_number(style_negative = "minus")) +
#     facet_grid(~ Treatment, space = "free", scales = "free") +
#     facetted_pos_scales(
#       x = list(
#         Treatment %in% c("Dark 15°C", "Light 15°C") ~
#           scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)),
#         Treatment %in% c("Light 20°C", "Light 25°C") ~
#           scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 30))
#       )
#     ) +
#     labs(x = "Detrital age (days)", y = "dm/dt") +
#     coord_cartesian(expand = F, clip = "off") +
#     mytheme
# 
# Fig_1c

require(ggridges)
# Fig_1c <- deco_prior_posterior_joined %>%
#   filter(Treatment != "Prior") %>%
#   select(starts_with("."), Treatment, k, tau) %>%
#   pivot_longer(cols = c(k, tau), names_to = "Parameter", values_to = "k") %>%
#   ggplot() +
#     stat_density_ridges(aes(k * 100, y = Parameter, fill = Treatment),
#                         from = 0, to = 16, n = 2^10, bandwidth = 15 * 0.02, 
#                         scale = 2, alpha = 0.8, colour = NA) +
#     scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
#                       guide = "none") +
#     facet_grid(~ Treatment) +
#     labs(x = expression("Decomposition (% d"^-1*")")) +
#     coord_cartesian(xlim = c(0, 15), expand = F, clip = "off") +
#     mytheme +
#     theme(axis.title.y = element_blank(),
#           axis.line.y = element_blank(),
#           axis.ticks.y = element_blank(),
#           axis.text.y = element_text(size = 12, vjust = -0.2, hjust = 0,
#                                      margin = margin(r = -0.64, unit = "cm")))
# 
# Fig_1c

Fig_1c <- deco_prior_posterior_joined %>%
  filter(Treatment != "Prior") %>%
  select(starts_with("."), Treatment, t0.5_k, t0.5_tau) %>%
  pivot_longer(cols = c(t0.5_k, t0.5_tau), names_to = "Parameter", 
               values_to = "t0.5", names_prefix = "t0.5_") %>%
  mutate(Parameter = Parameter %>% fct_relevel("tau")) %>%
  ggplot() +
    stat_density_ridges(aes(t0.5, y = Parameter, fill = Treatment),
                        from = 0, to = c(120, 120, 60, 60), n = 2^10, 
                        bandwidth = c(120, 120, 60, 60) * 0.02, 
                        scale = 1.5, alpha = 0.8, colour = NA) +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                      guide = "none") +
    scale_y_discrete(
      labels = c("k" = expression(italic("k")), "tau" = expression(italic("τ")))
    ) +
    facet_grid(~ Treatment, space = "free", scales = "free") +
    facetted_pos_scales(
      x = list(
        Treatment %in% c("Dark 15°C", "Light 15°C") ~
          scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)),
        Treatment %in% c("Light 20°C", "Light 25°C") ~
          scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 30))
      )
    ) +
    labs(x = "Half-life (days)") +
    coord_cartesian(expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_text(size = 12, vjust = -0.2, hjust = 1,
                                     margin = margin(r = 0.35, unit = "cm")))

Fig_1c


Fig_1 <- (
  ( Fig_1a + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.margin = margin(0, 0.5, 0.2, 0.2, unit = "cm")) ) / 
    ( Fig_1b + 
        theme(strip.text = element_blank(),
              plot.margin = margin(0.5, 0.5, 0.2, 0.2, unit = "cm")) ) / 
    ( Fig_1c + 
        theme(strip.text = element_blank()) )
) +
  plot_layout(heights = c(1, 1, 0.4))
  # plot_annotation(tag_levels = c("a", "b", "c")) &
  # theme(plot.tag = element_text(family = "Futura", size = 15, face = "bold"))

Fig_1

Fig_1 %>%
  ggsave(filename = "Fig_1.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 20, units = "cm")
