###################################
#### microdeco: photosynthesis ####
#### Luka Seamus Wright        ####
###################################

# 1. Prepare data ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)
oxygen <- here("Oxygen", "RDS", "Oxygen.rds") %>% 
  read_rds() %T>%
  print()

oxygen_summary <- oxygen %>%
  summarise(
    across(
      NP:S,
      list(mean = mean, sd = sd)
    ),
    # calculate probabilities that photosynthesis > 0
    NP_P = mean( NP > 0 ),
    GP_P = mean( GP > 0 ),
    n = n(),
    .by = ID:Mass
  ) %>%
  mutate(
    Tank = Tank %>% factor(),
    Treatment = case_when(
      Treatment == 1 ~ "Dark 15°C",
      Treatment == 2 ~ "Light 15°C",
      Treatment == 3 ~ "Light 20°C",
      Treatment == 4 ~ "Light 25°C"
    ) %>% fct_relevel("Dark 15°C")
  ) %>%
  select(!where(~ all(.x == 0))) %T>%
  print()

rm(oxygen)

deco <- here("Decomposition", "Decomposition.csv") %>% 
  read_csv() %>%
  mutate(
    Date = Date %>% dmy(),
    Day = Date[1] %--% Date / ddays(),
    Tank = Tank %>% factor(),
    Treatment = case_when(
      PAR == 0 ~ "Dark 15°C",
      Temperature == 15 ~ "Light 15°C",
      Temperature == 20 ~ "Light 20°C",
      Temperature == 25 ~ "Light 25°C"
    ) %>% fct_relevel("Dark 15°C")
  ) %>% 
  filter(!Day %in% c(0, 17, 20, 55)) %>% 
  # on days 17, 20 and 55, oxygen evolution was not measured
  select(Day, Tank, Treatment, Dry) %T>%
  print()

# 1.2 Calculate probability ####
prob <- oxygen_summary %>%
  select(Day, Tank, Treatment, NP_P, GP_P) %>%
  full_join(deco) %>%
  mutate(
    # replace NAs with 0s because when oxygen evolution was
    # not measured the tissue was either decomposed or not
    # intact, so can be assumed to have had no photosynthesis
    NP_P = if_else(is.na(NP_P), 0, NP_P),
    GP_P = if_else(is.na(GP_P), 0, GP_P),
    # reformat probabilities to have no 0s and 1s because the
    # beta distribution is not defined for these values
    NP_P = case_when(
      NP_P == 1 ~ 0.9999,
      NP_P == 0 ~ 0.0001,
      TRUE ~ NP_P
    ),
    GP_P = case_when(
      GP_P == 1 ~ 0.9999,
      GP_P == 0 ~ 0.0001,
      TRUE ~ GP_P
    )
  ) %>%
  arrange(Day, Tank) %T>%
  print()

rm(deco, oxygen_summary)

prob %>%
  write_rds(here("Oxygen", "RDS", "prob.rds"))

# 1.3 Explore data ####
# Define custom theme
mytheme <- theme(
  panel.background = element_blank(),
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
  text = element_text(family = "Futura")
)

prob %>%
  ggplot() +
    geom_point(
      aes(Day, NP_P), 
      shape = 16, size = 2, alpha = 0.4
    ) +
    facet_grid(~Treatment) +
    mytheme
# cf. Figure 4a in Wright et al. (2024), https://doi.org/10.1093/aob/mcad167

prob %>%
  ggplot() +
    geom_point(
      aes(Day, GP_P), 
      shape = 16, size = 2, alpha = 0.4
    ) +
    facet_grid(~Treatment) +
    mytheme
# cf. Figure 4b in Wright et al. (2024), https://doi.org/10.1093/aob/mcad167
# The probability of gross photosynthesis is more interesting because 
# no gross photosynthesis means that detritus is truly dead, so I will
# only analyse that.

# 2. Prior simulation ####
# I am using the model proposed by Wright (2026) QPB for detrital photosynthesis:
# p(t) = (alpha + tau) / ( 1 + exp( 5 / mu * (t - mu) ) ) - tau
# Since alpha = 1 and tau = 0 in this case where p is a probability, the model 
# simplifies to
# p(t) = 1 / ( 1 + exp( 5  /mu * (t - mu) ) ) = plogis( -5 / mu * (t - mu) )
# In contrast to the model used by Wright et al. (2024), this model fixes
# the intercept at plogis(5) = 0.99 because we know that detritus was alive at
# the start of the experiment. Hence, there is only one parameter: mu. This 
# is called the photosynthetic half-life. I am centring the prior on 50 as
# I have done in the decomposition model.

require(extraDistr)
tibble(
  n = 1:1e3,
  log_mu_mu = rnorm( 1e3 , log(50) , 0.3 ),
  # There are two factors in the multilevel model (treatment and tank)
  log_mu_sigma_t = rtnorm( 1e3 , 0 , 0.3 , 0 ),
  log_mu_sigma_ta = rtnorm( 1e3 , 0 , 0.3 , 0 ),
  nu = rgamma( 1e3 , 30^2/20^2 , 30/20^2 )
) %>%
  mutate(
    mu = exp( rnorm( n() , log_mu_mu , log_mu_sigma_t ) +
                rnorm( n() , 0 , log_mu_sigma_ta ) )
  ) %>%
  expand_grid(Day = prob %$% 
                seq(min(Day), max(Day), length.out = 100)) %>%
  mutate(
    P_mu = plogis( -5 / mu * ( Day - mu ) ),
    GP_P = rbeta( n() , P_mu * nu , (1 - P_mu) * nu )
  ) %>%
  pivot_longer(cols = c(P_mu, GP_P),
               names_to = "parameter") %>%
  ggplot(aes(Day, value, group = n)) +
    geom_line(alpha = 0.05) +
    coord_cartesian(expand = F, clip = "off") +
    facet_wrap(~parameter, scale = "free", nrow = 1) +
    theme_minimal() +
    theme(panel.grid = element_blank())

# 3. Stan model ####
require(cmdstanr)
prob_c_model <- here("Oxygen", "Stan", "prob_c.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

prob_nc_model <- here("Oxygen", "Stan", "prob_nc.stan") %>% 
  read_file() %>%
  write_stan_file() %>%
  cmdstan_model()

require(tidybayes)
prob_c_samples <- prob_c_model$sample(
          data = prob %>%
            select(Day, GP_P, Treatment, Tank) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

prob_nc_samples <- prob_nc_model$sample(
          data = prob %>%
            select(Day, GP_P, Treatment, Tank) %>%
            compose_data(),
          chains = 8,
          parallel_chains = parallel::detectCores(),
          iter_warmup = 1e4,
          iter_sampling = 1e4
        ) %T>%
  print()

# Save draws
prob_c_samples$draws() %>%
  write_rds(here("Oxygen", "RDS", "prob_c_samples.rds"))
prob_c_samples$draws(format = "df") %>%
  write_rds(here("Oxygen", "RDS", "prob_c_samples_df.rds"))

prob_nc_samples$draws() %>%
  write_rds(here("Oxygen", "RDS", "prob_nc_samples.rds"))
prob_nc_samples$draws(format = "df") %>%
  write_rds(here("Oxygen", "RDS", "prob_nc_samples_df.rds"))

# 4. Model checks ####
# 4.1 Rhat ####
prob_c_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# 8% of rhat above 1.001. rhat = 1.00 ± 0.000662.

prob_nc_samples$summary() %>%
  mutate(rhat_check = rhat > 1.001) %>%
  summarise(rhat_1.001 = sum(rhat_check) / length(rhat),
            rhat_mean = mean(rhat),
            rhat_sd = sd(rhat))
# No rhat above 1.001. rhat = 1.00 ± 0.0000755.
# Both models are fine but non-centred is the clear winner.

# 4.2 Chains ####
require(bayesplot)
prob_c_chains <- prob_c_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off") +
  mytheme

prob_c_chains %>%
  ggsave(filename = "prob_c_chains.pdf", path = here("Oxygen", "Plots"),
         device = cairo_pdf, width = 30, height = 20, units = "cm")

prob_nc_chains <- prob_nc_samples$draws(format = "df") %>%
  mcmc_rank_overlay() +
  guides(colour = guide_legend(nrow = 1)) +
  labs(title = "Non-centred model",
       y = "Frequency") +
  coord_cartesian(xlim = c(0, 8e4), ylim = c(0, 1e3),
                  expand = FALSE, clip = "off") +
  mytheme

prob_nc_chains %>%
  ggsave(filename = "prob_nc_chains.pdf", path = here("Oxygen", "Plots"),
         device = cairo_pdf, width = 30 * 7/5, height = 20 * 7/5, units = "cm")
# Chains are better for non-centred model

rm(prob_c_chains, prob_nc_chains) # Clean up
gc()

# 4.3 Pairs ####
prob_c_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("log_mu_mu", "log_mu_sigma_t", "log_mu_t[1]", "log_mu_t[2]",
             "log_mu_t[3]", "log_mu_t[4]", "log_mu_sigma_ta", "log_mu_ta[1]", 
             "log_mu_ta[2]", "log_mu_ta[5]", "log_mu_ta[6]", "log_mu_ta[10]", 
             "log_mu_ta[11]", "log_mu_ta[13]", "log_mu_ta[14]", "nu"),
    grid_args = list(top = "Centred model")
  ) %>%
  ggsave(filename = "prob_c_pairs.png", path = here("Oxygen", "Plots"),
         width = 100, height = 100, units = "cm", bg = "white")

prob_nc_samples$draws(format = "df") %>%
  mcmc_pairs(
    pars = c("log_mu_mu", "log_mu_sigma_t", "log_mu_t[1]", "log_mu_t[2]",
             "log_mu_t[3]", "log_mu_t[4]", "log_mu_sigma_ta", "log_mu_ta[1]", 
             "log_mu_ta[2]", "log_mu_ta[5]", "log_mu_ta[6]", "log_mu_ta[10]", 
             "log_mu_ta[11]", "log_mu_ta[13]", "log_mu_ta[14]", "nu"),
    grid_args = list(top = "Non-centred model")
  ) %>%
  ggsave(filename = "prob_nc_pairs.png", path = here("Oxygen", "Plots"),
         width = 100, height = 100, units = "cm", bg = "white")
# Pairs look fine for both models.

# 5. Prior-posterior comparison ####
# Hierarchical priors cannot effectively sampled when centred.
# Hence the non-centred model will be used to sample priors.
source("functions.R")
prob_prior <- prior_samples(
  model = prob_nc_model,
  data = prob %>%
    select(Day, GP_P, Treatment, Tank) %>%
    compose_data()
)

# Centred model
prob_c_prior_posterior_treatment <- prob_prior %>% 
  prior_posterior_draws(
    posterior_samples = prob_c_samples,
    group = prob %>% select(Treatment),
    parameters = c("log_mu_mu", "log_mu_sigma_t", 
                   "log_mu_t[Treatment]",
                   "log_mu_sigma_ta", "nu"),
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

prob_c_prior_posterior_tank <- prob_prior %>% 
  prior_posterior_draws(
    posterior_samples = prob_c_samples,
    group = prob %>% select(Tank),
    parameters = c("log_mu_ta[Tank]"),
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
prob_nc_prior_posterior_treatment <- prob_prior %>% 
  prior_posterior_draws(
    posterior_samples = prob_nc_samples,
    group = prob %>% select(Treatment),
    parameters = c("log_mu_mu", "log_mu_sigma_t", 
                   "log_mu_t[Treatment]",
                   "log_mu_sigma_ta", "nu"),
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

prob_nc_prior_posterior_tank <- prob_prior %>% 
  prior_posterior_draws(
    posterior_samples = prob_nc_samples,
    group = prob %>% select(Tank),
    parameters = c("log_mu_ta[Tank]"),
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
prob_prior_posterior_plot <- 
  ( prob_c_prior_posterior_treatment / prob_c_prior_posterior_tank +
        plot_layout(heights = c(1, 0.8)) ) | 
  ( prob_nc_prior_posterior_treatment / prob_nc_prior_posterior_tank ) +
        plot_layout(heights = c(1, 0.8))

prob_prior_posterior_plot %>%
  ggsave(filename = "prob_prior_posterior.pdf", 
         path = here("Oxygen", "Plots"),
         device = cairo_pdf, width = 50, height = 30, units = "cm")
# Model posteriors look nearly identical. Choose non-centred as optimal model.

rm(prob_nc_prior_posterior_treatment, prob_nc_prior_posterior_tank,
   prob_c_prior_posterior_treatment, prob_c_prior_posterior_tank,
   prob_prior_posterior_plot, prob_c_model, prob_c_samples) # Clean up
gc()

# 6. Photosynthetic half-life ####
# Global (across treatments and tanks)
prob_prior_posterior_global <- prob_prior %>% 
  prior_posterior_draws(
    posterior_samples = prob_nc_samples,
    parameters = c("log_mu_mu", "log_mu_sigma_t", 
                   "log_mu_sigma_ta", "nu"),
    format = "short"
  ) %>%
  mutate( # Calculate mu for unobserved treatments and tanks
    mu = exp(
      rnorm( n() , log_mu_mu , log_mu_sigma_t ) +
        rnorm( n() , 0 , log_mu_sigma_ta )
    )
  ) %>%
  select(starts_with("."), distribution, mu, nu) %T>%
  print()

prob_prior_posterior_global %>%
  pivot_longer(cols = -c(starts_with("."), distribution),
               names_to = "parameter") %>%
  summarise(mean = mean(value), sd = sd(value), n = n(),
            .by = c(distribution, parameter)) %>%
  arrange(parameter)

# Treatments (across tanks)
prob_prior_posterior <- prob_prior %>% 
  prior_posterior_draws(
    posterior_samples = prob_nc_samples,
    group = prob %>% select(Treatment),
    parameters = c("log_mu_t[Treatment]", "log_mu_sigma_ta", "nu"),
    format = "short"
  ) %>% 
  mutate( # Calculate mu for new, unobserved tanks
    mu = exp( rnorm( n() , log_mu_t , log_mu_sigma_ta ) )
  ) %>% # Remove redundant priors (keep only dark prior)
  filter(!(Treatment %>% str_detect("Light") & distribution == "prior")) %>%
  mutate( # Embed prior in treatment
    Treatment = if_else(
      distribution == "prior", "Prior", Treatment
    ) %>% fct()
  ) %>%
  select(starts_with("."), Treatment, mu, nu) %T>%
  print()

prob_prior_posterior %>%
  pivot_longer(cols = -c(starts_with("."), Treatment),
               names_to = "parameter") %>%
  summarise(mean = mean(value), sd = sd(value), n = n(),
            .by = c(Treatment, parameter)) %>%
  arrange(parameter)

# Save parameter distributions
prob_prior_posterior_global %>%
  write_rds(here("Oxygen", "RDS", "prob_prior_posterior_global.rds"))
prob_prior_posterior %>%
  write_rds(here("Oxygen", "RDS", "prob_prior_posterior.rds"))

# 7. Prediction ####
# Predict across predictor range
prob_prediction <- prob_prior_posterior %>%
  spread_continuous(
    data = prob,
    group_name = "Treatment",
    predictor_name = "Day"
  ) %>%
  mutate(
    P_mu = plogis( -5 / mu * ( Day - mu ) ),
    GP_P = rbeta( n() , P_mu * nu , (1 - P_mu) * nu )
  ) %T>%
  print()

# Summarise
prob_prediction_summary <- prob_prediction %>%
  summarise(
    across(
      c(P_mu, GP_P),
      list(mean = mean,
           median = median, 
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
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

# Clean up
rm(prob_prediction)
gc()

# Save predictions
prob_prediction_summary %>%
  write_rds(here("Oxygen", "RDS", "prob_prediction.rds"))

# 8. Tables ####
# 8.1 Load posteriors for detrital half-life ####
deco_prior_posterior <- here("Decomposition", "RDS", "deco_prior_posterior.rds") %>%
  read_rds() %T>%
  print()

deco_k_prior_posterior <- here("Decomposition", "RDS", "deco_k_prior_posterior.rds") %>%
  read_rds() %T>%
  print()

# 8.2 Join posteriors ####
prior_posterior_joined <- prob_prior_posterior %>% 
  select(starts_with("."), Treatment, mu) %>%
  full_join(
    deco_prior_posterior %>% 
      select(starts_with("."), Treatment, tau)
  ) %>%
  full_join(
    deco_k_prior_posterior %>% 
      select(starts_with("."), Treatment, k)
  ) %>%
  # Calculate detrital half-lives based on k and tau (i.e. final k),
  # as well as differences between detrital and photosynthetic half-lives
  mutate(
    t0.5_k = log(2)/k,
    t0.5_tau = log(2)/tau,
    k_tau_diff = t0.5_k - t0.5_tau,
    mu_k_diff = mu - t0.5_k
  ) %T>%
  print()

# 8.3 Summarise half-lives ####
require(glue)
Table_1 <- prior_posterior_joined %>%
  mutate(
    k = k * 100, # Convert exponential rates to %
    tau = tau * 100
  ) %>%
  summarise(
    across(
      c(k, tau, mu, t0.5_k, t0.5_tau, k_tau_diff, mu_k_diff),
      list(mean = mean, sd = sd, median = median)
    ),
    P_k_tau = mean( k_tau_diff > 0 ),
    P_mu_k = mean( mu_k_diff > 0 ),
    n = n(),
    .by = Treatment
  ) %>%
  mutate(
    across(where(is.numeric), ~signif(.x, 2)),
    k = glue("{k_mean} ± {k_sd} ({k_median})"),
    tau = glue("{tau_mean} ± {tau_sd} ({tau_median})"),
    mu = glue("{mu_mean} ± {mu_sd} ({mu_median})"),
    t0.5_k = glue("{t0.5_k_mean} ± {t0.5_k_sd} ({t0.5_k_median})"),
    t0.5_tau = glue("{t0.5_tau_mean} ± {t0.5_tau_sd} ({t0.5_tau_median})"),
    k_tau_diff = glue("{k_tau_diff_mean} ± {k_tau_diff_sd} ({k_tau_diff_median})"),
    mu_k_diff = glue("{mu_k_diff_mean} ± {mu_k_diff_sd} ({mu_k_diff_median})")
  ) %>%
  select(!c(ends_with("mean"), ends_with("median"), ends_with("sd"))) %T>%
  print()

Table_1 %>%
  write_csv(here("Tables", "Table_1.csv"))

require(officer)
read_docx() %>%
  body_add_table(value = Table_1) %>%
  print(target = here("Tables", "Table_1.docx"))

# 8.4 Treatment contrasts ####
prob_contrast <- prob_prior_posterior %>%
  filter(Treatment != "Prior") %>%
  droplevels() %>%
  select(starts_with("."), Treatment, mu) %>%
  pivot_wider(names_from = Treatment,
              values_from = mu) %>%
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

prob_contrast_summary <- prob_contrast %>%
  summarise(
    across(
      c(diff, ratio),
      list(mean = mean, sd = sd, median = median)
    ),
    P = max( mean( diff > 0 ) , mean( diff < 0 ) ),
    n = n(),
    .by = contrast
  ) %>%
  mutate(
    across(where(is.numeric), ~signif(.x, 2)),
    diff = glue("{diff_mean} ± {diff_sd} ({diff_median})"),
    ratio = glue("{ratio_mean} ± {ratio_sd} ({ratio_median})")
  ) %>%
  select(!c(ends_with("mean"), ends_with("median"), ends_with("sd"))) %T>%
  print()

prob_contrast_summary %>%
  write_csv(here("Tables", "prob_contrast.csv"))

read_docx() %>%
  body_add_table(value = prob_contrast_summary) %>%
  print(target = here("Tables", "prob_contrast.docx"))

# 9. Figure ####
# 9.1 Load predictions for decomposition ####
deco_prediction <- here("Decomposition", "RDS", "deco_prediction.rds") %>%
  read_rds() %T>%
  print()

deco_summary <- here("Decomposition", "RDS", "deco_summary.rds") %>%
  read_rds() %T>%
  print()

# 9.2 First panel ####
require(ggh4x)
Fig_1a <- prob_prediction_summary %>%
  filter(Treatment != "Prior") %>%
  ggplot() +
    geom_point(
      data = prob,
      aes(Day, GP_P, colour = Treatment),
      shape = 16, alpha = 0.5, size = 2.35 # this matches default size of point range
    ) +
    geom_line(aes(Day, P_mu.median, colour = Treatment)) +
    geom_ribbon(aes(Day, ymin = P_mu.lower, ymax = P_mu.upper,
                    alpha = factor(.width), fill = Treatment)) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                        guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                      guide = "none") +
    scale_alpha_manual(values = c(0.5, 0.4, 0.3), guide = "none") +
    scale_y_continuous(labels = scales::label_number(accuracy = c(1, 1e-2, 0.1, 1e-2, 1))) +
    facet_grid(~ Treatment, space = "free", scales = "free") +
    facetted_pos_scales(
      x = list(
        Treatment %in% c("Dark 15°C", "Light 15°C") ~
          scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)),
        Treatment %in% c("Light 20°C", "Light 25°C") ~
          scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 30))
      )
    ) +
    labs(x = "Detrital age (days)", 
         y = expression("Photosynthesis ("*italic("P")*")")) +
    coord_cartesian(ylim = c(0, 1), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_text(margin = margin(r = -0.3, unit = "cm")))

Fig_1a

# 9.3 Second panel ####
Fig_1b <- deco_prediction %>%
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
    labs(x = "Detrital age (days)", y = "Detrital mass (%)") +
    coord_cartesian(ylim = c(0, 160), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_text(margin = margin(r = -0.3, unit = "cm")))

Fig_1b

# 9.4 Third panel ####
Fig_1c <- deco_prediction %>%
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

Fig_1c

# 9.5 Fourth panel ####
require(ggridges)
Fig_1d <- prior_posterior_joined %>%
  filter(Treatment != "Prior") %>%
  select(starts_with("."), Treatment, mu, t0.5_k, t0.5_tau) %>%
  pivot_longer(cols = c(mu, t0.5_k, t0.5_tau), names_to = "Parameter", 
               values_to = "t0.5", names_prefix = "t0.5_") %>%
  mutate(Parameter = Parameter %>% fct_relevel("tau")) %>%
  ggplot() +
    stat_density_ridges(aes(t0.5, y = Parameter, fill = Treatment),
                        from = 0, to = c(140, 140, 80, 80), n = 2^10, 
                        bandwidth = c(140, 140, 80, 80) * 0.02, 
                        scale = 1.5, alpha = 0.8, colour = NA) +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                      guide = "none") +
    scale_y_discrete(
      labels = c("mu" = expression(italic("μ")), 
                 "k" = expression(italic("t")["½ "[italic("k")]]), 
                 "tau" = expression(italic("t")["½ "[italic("τ")]]))
      # Greek letters show as bold in Futura, so this must be changed later.
      # Sub-subscript k and tau must later be adjusted to the same size as ½.
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
          axis.text.y = element_text(size = 12, vjust = 0, hjust = 0))

Fig_1d

# 9.6 Combined figure ####
Fig_1 <- (
  ( Fig_1a + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.margin = margin(0, 0.5, 0.2, 0.2, unit = "cm")) ) / 
    ( Fig_1b + 
        theme(strip.text = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              plot.margin = margin(0.5, 0.5, 0.2, 0.2, unit = "cm")) ) / 
    ( Fig_1c + 
        theme(strip.text = element_blank(),
              plot.margin = margin(0.5, 0.5, 0.2, 0.2, unit = "cm")) ) / 
    ( Fig_1d + 
        theme(strip.text = element_blank()) )
) +
  plot_layout(heights = c(1, 1, 1, 0.8))

Fig_1

Fig_1 %>%
  ggsave(filename = "Fig_1.pdf", path = "Figures",
         device = cairo_pdf, height = 20, width = 20, units = "cm")