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

# 1.3 Explore data ####
prob %>%
  ggplot() +
    geom_point(
      aes(Day, NP_P), 
      shape = 16, size = 2, alpha = 0.4
    ) +
    facet_grid(~Treatment) +
    theme_minimal()
# cf. Figure 4a in Wright et al. (2024), https://doi.org/10.1093/aob/mcad167

prob %>%
  ggplot() +
    geom_point(
      aes(Day, GP_P), 
      shape = 16, size = 2, alpha = 0.4
    ) +
    facet_grid(~Treatment) +
    theme_minimal()
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


# 5. Prior-posterior comparison ####

# 6. Photosynthetic half-life ####

# 7. Prediction ####

# 8. Figure ####

