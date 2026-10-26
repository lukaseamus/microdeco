# 1. Photosynthesis ####
# 1.1 Load data ####
require(tidyverse)
require(magrittr)
require(here)

prob <- here("Oxygen", "RDS", "prob.rds") %>%
  read_rds() %T>%
  print()

prob_prior_posterior <- here("Oxygen", "RDS", "prob_prior_posterior.rds") %>%
  read_rds() %T>%
  print()

# Predict across full predictor range
source("functions.R")
require(tidybayes)
prob_prediction <- prob_prior_posterior %>%
  spread_continuous(
    data = prob,
    predictor_name = "Day",
    length = 600 # increase resolution because prior has a huge range
  ) %>%
  mutate(
    P_mu = plogis( -5 / mu * ( Day - mu ) ),
    GP_P = rbeta( n() , P_mu * nu , (1 - P_mu) * nu )
  ) %>%
  # below is a longer but much much (!) more efficient version of 
  # group_by(Day, Treatment) %>% 
  # tidybayes::median_qi(P_mu, GP_P, .width = c(.5, .8, .9))
  summarise(
    across(
      c(P_mu, GP_P),
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
  rename(P_mu = P_mu.median, GP_P = GP_P.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

# 1.2 Aesthetics ####
# Data points
prob_points <- prob %>%
  mutate(
    alpha = 0.5, # Baseline transparency
    colour = case_when( # Colours
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %>%
  rename(x = Day, y = GP_P) %>% # Rename to x and y for consistency
  select(-c(Tank, NP_P, Dry)) %T>% # Remove irrelevant variables
  print()

# Lines
prob_line <- prob_prediction %>%
  select(Treatment, Day, P_mu) %>%
  rename(x = Day, y = P_mu) %>%
  mutate(
    colour = case_when(
      Treatment == "Prior" ~ "#b5b8ba", # here I also have Prior
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %T>%
  print()

# Ribbons
prob_ribbon <- prob_prediction %>%
  group_by(Treatment, .width) %>%
  reframe(
    x = c(Day, rev(Day)),
    y = c(P_mu.upper, rev(P_mu.lower))
  ) %>%
  mutate(
    alpha = case_when(
      .width == .9 ~ .5, 
      .width == .8 ~ .4, 
      .width == .5 ~ .3
    ),
    fill = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %>%
  select(-.width) %T>%
  print()

# Densities
prob_dens <- prob_prior_posterior %>%
  group_by(Treatment) %>% 
  reframe( 
    # Close polygons at zero and max with c(0, density, max) for x and
    # at zero and zero with c(0, density, 0) for y.
    x = c(0, density(mu, n = 2^10, from = 0, to = 150, bw = 150 * 0.01)$x, 150),
    y = c(0, density(mu, n = 2^10, from = 0, to = 150, bw = 150 * 0.01)$y, 0)
  ) %>%
  group_by(Treatment) %>% # Standardise area with Riemann sum (avoid manually added x[1]).
  mutate( y = y * 3 / ( sum(y) * ( x[3] - x[2] ) ) ) %>%
  ungroup() %T>%
  print()

prob_dens %<>%
  mutate(
    fill = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %T>%
  print()

# 1.3 Static figures ####
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

# Top panel
ggplot() +
  geom_point(
    data = prob_points,
    aes(x, y, colour = colour, alpha = alpha),
    shape = 16, size = 2.35 # this size matches default size of point range
  ) +
  geom_line(
    data = prob_line,
    aes(x, y, colour = colour)
  ) +
  geom_polygon(
    data = prob_ribbon,
    aes(x, y, fill = fill, alpha = alpha, group = alpha)
  ) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::label_number(accuracy = c(1, 1e-2, 0.1, 1e-2, 1))) +
  facet_grid(~ Treatment, space = "free", scales = "free") +
  labs(x = "Detrital age (days)", 
       y = expression("Photosynthesis ("*italic("P")*")")) +
  coord_cartesian(xlim = c(0, 120), ylim = c(0, 1), expand = F, clip = "off") +
  mytheme

# Bottom panel
ggplot() +
  geom_polygon(
    data = prob_dens,
    aes(x = x, y = y, fill = fill)
  ) +
  scale_fill_identity() +
  # Ensure same spacing as top panel
  scale_y_continuous(labels = scales::label_number(accuracy = 1e-2)) +
  facet_grid(~ Treatment, space = "free", scales = "free") +
  labs(x = "Half-life (days)",
       y = expression("Photosynthesis ("*italic("P")*")")) +
  coord_cartesian(xlim = c(0, 120), expand = F, clip = "off") +
  mytheme +
  # These will be hidden (transparent) but ensure identical spacing
  theme(axis.title.y = element_text(colour = "transparent"),
        axis.text.y = element_text(colour = "transparent"),
        axis.ticks.y = element_line(colour = "transparent"),
        axis.line.y = element_line(colour = "transparent"))

# 1.4 Animation ####
# Define enter/exit functions
enter <- function(data){
  data$alpha <- 0
  data
}

exit <- function(data){
  data$alpha <- 0
  data
}

# Tween points
require(tweenr)
prob_points_ani <- bind_rows( 
  tween_state( # Points are arbitrarily paired
    prob_points %>% filter(Treatment == "Prior"), 
    prob_points %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    prob_points %>% filter(Treatment == "Dark 15°C"), 
    prob_points %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_state(
    prob_points %>% filter(Treatment == "Light 15°C"), 
    prob_points %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_state(
    prob_points %>% filter(Treatment == "Light 20°C"), 
    prob_points %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_state(
    prob_points %>% filter(Treatment == "Light 25°C"), 
    prob_points %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_state(
    prob_points %>% filter(Treatment == "Light 15°C"), 
    prob_points %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_state(
    prob_points %>% filter(Treatment == "Dark 15°C"), 
    prob_points %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_state(
    prob_points %>% filter(Treatment == "Light 20°C"), 
    prob_points %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150),
  tween_state(
    prob_points %>% filter(Treatment == "Light 25°C"), 
    prob_points %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 8 * 150)
) %T>%
  print()

# Tween lines
# tween_state also works for some lines, but 
# I believe it is more stable to use tween_path.
require(transformr)
prob_line_ani <- bind_rows(
  tween_path(
    prob_line %>% filter(Treatment == "Prior"), 
    prob_line %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_path(
    prob_line %>% filter(Treatment == "Dark 15°C"), 
    prob_line %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_path(
    prob_line %>% filter(Treatment == "Light 15°C"), 
    prob_line %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_path(
    prob_line %>% filter(Treatment == "Light 20°C"), 
    prob_line %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_path(
    prob_line %>% filter(Treatment == "Light 25°C"), 
    prob_line %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_path(
    prob_line %>% filter(Treatment == "Light 15°C"), 
    prob_line %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_path(
    prob_line %>% filter(Treatment == "Dark 15°C"), 
    prob_line %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_path(
    prob_line %>% filter(Treatment == "Light 20°C"), 
    prob_line %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150),
  tween_path(
    prob_line %>% filter(Treatment == "Light 25°C"), 
    prob_line %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 8 * 150)
) %T>%
  print()

# Tween ribbons
# tween_path is more stable than tween_polygon for ribbons.
prob_ribbon_ani <- bind_rows(
  tween_path(
    prob_ribbon %>% filter(Treatment == "Prior"),
    prob_ribbon %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_path(
    prob_ribbon %>% filter(Treatment == "Dark 15°C"),
    prob_ribbon %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_path(
    prob_ribbon %>% filter(Treatment == "Light 15°C"),
    prob_ribbon %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_path(
    prob_ribbon %>% filter(Treatment == "Light 20°C"),
    prob_ribbon %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_path(
    prob_ribbon %>% filter(Treatment == "Light 25°C"),
    prob_ribbon %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_path(
    prob_ribbon %>% filter(Treatment == "Light 15°C"),
    prob_ribbon %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_path(
    prob_ribbon %>% filter(Treatment == "Dark 15°C"),
    prob_ribbon %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_path(
    prob_ribbon %>% filter(Treatment == "Light 20°C"),
    prob_ribbon %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150),
  tween_path(
    prob_ribbon %>% filter(Treatment == "Light 25°C"),
    prob_ribbon %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 8 * 150)
) %T>%
  print()

# Tween densities
prob_dens_ani <- bind_rows(
  tween_polygon(
    prob_dens %>% filter(Treatment == "Prior"),
    prob_dens %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Dark 15°C"),
    prob_dens %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 15°C"),
    prob_dens %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 20°C"),
    prob_dens %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 25°C"),
    prob_dens %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 15°C"),
    prob_dens %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Dark 15°C"),
    prob_dens %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 20°C"),
    prob_dens %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 25°C"),
    prob_dens %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 8 * 150)
) %T>%
  print()

# 1.5 Animated figures ####
require(gganimate) # also install gifski
require(ggtext) # for html syntax in gif

# Top panel
( ggplot() +
    geom_point(
      data = prob_points_ani,
      aes(x, y, colour = colour, alpha = alpha),
      shape = 16, size = 2.35 # this size matches default size of point range
    ) +
    geom_line(
      data = prob_line_ani,
      aes(x, y, colour = colour)
    ) +
    geom_polygon(
      data = prob_ribbon_ani,
      aes(x, y, fill = fill, alpha = alpha, group = alpha)
    ) +
    geom_text(
      data = prob_line_ani %>% distinct(Treatment, colour, .frame),
      aes(x = 0.8, y = 0.05, label = Treatment, colour = colour),
      hjust = 0, vjust = 0, size.unit = "pt", size = 12, family = "Futura"
    ) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_y_continuous(labels = scales::label_number(accuracy = c(1, 1e-2, 0.1, 1e-2, 1))) +
    labs(
      x = "Detrital age (days)", 
      y = "Photosynthesis (*P*)"
    ) +
    coord_cartesian(xlim = c(0, 120), ylim = c(0, 1), expand = F, clip = "off") +
    transition_manual(.frame) +
    mytheme +
    theme(
      axis.title.y = element_markdown() # only specify for y axis
    ) ) %>%
  animate(
    fps = 25, # optimal gif frame rate (4 centiseconds per frame) with reasonable file size
    duration = 9 * 3, # number of transitions (9) * time per transition (3 s) to set speed
    detail = 2, # 25*9*3 is 1/2 of all frames (9*150), so the skipped frames can add detail
    width = 20, height = 6,
    units = "cm", res = 300,
    renderer = gifski_renderer(),
    device = "ragg_png" # necessary for italic!
  ) %>%
  anim_save(
    filename = "prob_top.gif", 
    path = here("Figures", "Animation")
  )

# Bottom panel
( ggplot() +
    geom_polygon(
      data = prob_dens_ani,
      aes(x = x, y = y, fill = fill)
    ) +
    scale_fill_identity() +
    # Ensure same spacing as top panel
    scale_y_continuous(labels = scales::label_number(accuracy = 1e-2)) +
    labs(
      x = "Photosynthetic half-life (days)",
      y = "Photosynthesis (*P*)"
    ) +
    coord_cartesian(xlim = c(0, 120), expand = F, clip = "off") +
    transition_manual(.frame) +
    mytheme +
    theme(
      # These will be hidden (transparent) but ensure identical spacing
      axis.title.y = element_markdown(colour = "transparent"),
      axis.text.y = element_text(colour = "transparent"),
      axis.ticks.y = element_line(colour = "transparent"),
      axis.line.y = element_line(colour = "transparent")
    ) ) %>%
  animate(
    fps = 25, duration = 9 * 3, detail = 2,
    width = 20, height = 4, units = "cm", res = 300, 
    renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "prob_bottom.gif", 
    path = here("Figures", "Animation")
  )

# 2. Decomposition ####
# 2.1 Load data ####
deco <- here("Decomposition", "RDS", "deco_summary.rds") %>%
  read_rds() %T>%
  print()

deco_prior_posterior <- here("Decomposition", "RDS", "deco_prior_posterior.rds") %>%
  read_rds() %T>%
  print()

deco_k_prior_posterior <- here("Decomposition", "RDS", "deco_k_prior_posterior.rds") %>%
  read_rds() %T>%
  print()

# R doesn't have a built-in log1p_exp function
log1p_exp <- function(x) {
  ifelse(
    x > 0, 
    x + log1p(exp(-x)),
    log1p(exp(x))
  )
}

# Predict across full predictor range
require(extraDistr) # needed for beta prime
deco_prediction <- deco_prior_posterior %>%
  spread_continuous(
    data = deco %>%
      add_row(
        Day = 0, # Ensure predictor range starts at 0
        Treatment = deco %$% levels(Treatment) %>% fct()
      ),
    predictor_name = "Day",
    length = 200 # double resolution
  ) %>%
  mutate(
    r_mu = exp(
      Day * alpha - ( alpha + tau ) * mu / 5 * (
        log1p_exp( 5 / mu * ( Day - mu ) ) - log1p_exp( -5 )
      )
    ),
    nu = ( epsilon - theta ) * exp( -lambda * Day ) + theta,
    r = rbetapr( n() , r_mu * ( 1 + nu ) , 2 + nu )
  ) %>%
  summarise(
    across(
      c(r_mu, r),
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
  rename(r_mu = r_mu.median, r = r.median) %>%
  pivot_longer(cols = contains("lower") | contains("upper")) %>%
  separate(col = name, into = c("name", ".width"), sep = "_(?=[^_]*$)") %>%
  pivot_wider(names_from = name, values_from = value) %T>%
  print()

# 2.2 Aesthetics ####
# Data points
deco_points <- deco %>%
  mutate(
    alpha = 0.5, # Baseline transparency
    colour = case_when( # Colours
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %>%
  rename(x = Day, y = Ratio_mean) %>% # Rename to x and y for consistency
  select(Treatment, x, y, alpha, colour) %T>% # Select relevant variables
  print()

# Error bars
deco_bars <- deco %>%
  mutate(Replicate = row_number()) %>%
  select(Treatment, Replicate, Day, Ratio_lower, Ratio_upper) %>%
  pivot_longer(cols = starts_with("Ratio"), values_to = "y") %>%
  rename(x = Day) %>%
  mutate(
    alpha = 0.5, # Baseline transparency
    colour = case_when( # Colours
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %>%
  select(-name) %T>%
  print()

# Lines
deco_line <- deco_prediction %>%
  select(Treatment, Day, r_mu) %>%
  rename(x = Day, y = r_mu) %>%
  mutate(
    colour = case_when(
      Treatment == "Prior" ~ "#b5b8ba", # here I also have Prior
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %T>%
  print()

# Ribbons
deco_ribbon <- deco_prediction %>%
  group_by(Treatment, .width) %>%
  reframe(
    x = c(Day, rev(Day)),
    y = c(r_mu.upper, rev(r_mu.lower))
  ) %>%
  mutate(
    alpha = case_when(
      .width == .9 ~ .5, 
      .width == .8 ~ .4, 
      .width == .5 ~ .3
    ),
    fill = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %>%
  select(-.width) %T>%
  print()

# Densities
deco_dens <- deco_k_prior_posterior %>%
  mutate(t0.5 = log(2)/k) %>% # Calculate detrital half-life (t½)
  group_by(Treatment) %>% 
  reframe(
    x = c(0, density(t0.5, n = 2^10, from = 0, to = 150, bw = 150 * 0.01)$x, 150),
    y = c(0, density(t0.5, n = 2^10, from = 0, to = 150, bw = 150 * 0.01)$y, 0)
  ) %>%
  group_by(Treatment) %>% # note I use same multiplier (3) as for prob_dens 
  mutate( y = y * 3 / ( sum(y) * ( x[3] - x[2] ) ) ) %>%
  ungroup() %T>%
  print()

deco_dens %<>%
  mutate(
    fill = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %T>%
  print()

# 2.3 Static figures ####
# Top panel
ggplot() +
  geom_hline(yintercept = 1) +
  geom_point(
    data = deco_points,
    aes(x, y, colour = colour, alpha = alpha),
    shape = 16, size = 2.35
  ) +
  geom_line(
    data = deco_bars,
    aes(x, y, colour = colour, 
        alpha = alpha, group = Replicate),
  ) +
  geom_line(
    data = deco_line,
    aes(x, y, colour = colour)
  ) +
  geom_polygon(
    data = deco_ribbon,
    aes(x, y, fill = fill, alpha = alpha, group = alpha)
  ) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  scale_y_continuous(labels = scales::label_number(accuracy = c(1, rep(0.1, 4)))) +
  facet_grid(~ Treatment, space = "free", scales = "free") +
  labs(x = "Detrital age (days)", 
       y = expression("Detrital mass ("*italic("m")*"/"*italic("m")[0]*")")) +
  coord_cartesian(xlim = c(0, 120), ylim = c(0, 1.6), expand = F, clip = "off") +
  mytheme

# Bottom panel v1
ggplot() +
  geom_polygon(
    data = deco_dens,
    aes(x = x, y = y, fill = fill)
  ) +
  scale_fill_identity() +
  # Ensure same spacing as top panel
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  facet_grid(~ Treatment, space = "free", scales = "free") +
  labs(x = "Half-life (days)",
       y = expression("Detrital mass ("*italic("m")*"/"*italic("m")[0]*")")) +
  coord_cartesian(xlim = c(0, 120), ylim = c(0, 0.25), expand = F, clip = "off") +
  mytheme +
  # These will be hidden (transparent) but ensure identical spacing
  theme(axis.title.y = element_text(colour = "transparent"),
        axis.text.y = element_text(colour = "transparent"),
        axis.ticks.y = element_line(colour = "transparent"),
        axis.line.y = element_line(colour = "transparent"))

# Bottom panel v2
ggplot() +
  geom_polygon(
    data = prob_dens,
    aes(x = x, y = y + 0.1, fill = fill), alpha = 0.8
  ) +
  geom_polygon(
    data = deco_dens,
    aes(x = x, y = y, fill = fill), alpha = 0.8
  ) +
  geom_text(
    data = tibble(
      Treatment = deco_dens %$% levels(Treatment) %>% 
        fct() %>% rep(each = 2),
      label = c("Photosynth.", "Detrital", rep(NA, 8)),
      x = -100, y = c(0.1, 0, rep(NA, 8))
    ),
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 0, family = "Futura",
    size.unit = "pt", size = 12
  ) +
  scale_fill_identity() +
  # Ensure same spacing as top panel
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  facet_grid(~ Treatment, space = "free", scales = "free") +
  labs(x = "Half-life (days)",
       y = expression("Detrital mass ("*italic("m")*"/"*italic("m")[0]*")")) +
  coord_cartesian(xlim = c(0, 120), expand = F, clip = "off") +
  mytheme +
  # These will be hidden (transparent) but ensure identical spacing
  theme(axis.title.y = element_text(colour = "transparent"),
        axis.text.y = element_text(colour = "transparent"),
        axis.ticks.y = element_line(colour = "transparent"),
        axis.line.y = element_line(colour = "transparent"))

# or even simpler:
ggplot() +
  geom_polygon(
    data = prob_dens,
    aes(x = x, y = y, fill = fill), alpha = 0.8
  ) +
  geom_polygon(
    data = deco_dens,
    aes(x = x, y = y, fill = fill), alpha = 0.8
  ) +
  scale_fill_identity() +
  # Ensure same spacing as top panel
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
  facet_grid(~ Treatment, space = "free", scales = "free") +
  labs(x = "Half-life (days)",
       y = expression("Detrital mass ("*italic("m")*"/"*italic("m")[0]*")")) +
  coord_cartesian(xlim = c(0, 120), ylim = c(0, 0.25), expand = F, clip = "off") +
  mytheme +
  # These will be hidden (transparent) but ensure identical spacing
  theme(axis.title.y = element_text(colour = "transparent"),
        axis.text.y = element_text(colour = "transparent"),
        axis.ticks.y = element_line(colour = "transparent"),
        axis.line.y = element_line(colour = "transparent"))

# 2.4 Animation ####
# Tween points
deco_points_ani <- bind_rows( 
  tween_state( # Points are arbitrarily paired
    deco_points %>% filter(Treatment == "Prior"), 
    deco_points %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    deco_points %>% filter(Treatment == "Light 15°C"), 
    deco_points %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_state(
    deco_points %>% filter(Treatment == "Dark 15°C"), 
    deco_points %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_state(
    deco_points %>% filter(Treatment == "Light 20°C"), 
    deco_points %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_state(
    deco_points %>% filter(Treatment == "Light 25°C"), 
    deco_points %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween bars
deco_bars_ani <- bind_rows( 
  tween_state( # Bars are treated the same as points
    deco_bars %>% filter(Treatment == "Prior"), 
    deco_bars %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    deco_bars %>% filter(Treatment == "Light 15°C"), 
    deco_bars %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_state(
    deco_bars %>% filter(Treatment == "Dark 15°C"), 
    deco_bars %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_state(
    deco_bars %>% filter(Treatment == "Light 20°C"), 
    deco_bars %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_state(
    deco_bars %>% filter(Treatment == "Light 25°C"), 
    deco_bars %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween lines
deco_line_ani <- bind_rows(
  tween_path(
    deco_line %>% filter(Treatment == "Prior"), 
    deco_line %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_path(
    deco_line %>% filter(Treatment == "Light 15°C"), 
    deco_line %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_path(
    deco_line %>% filter(Treatment == "Dark 15°C"), 
    deco_line %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_path(
    deco_line %>% filter(Treatment == "Light 20°C"), 
    deco_line %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_path(
    deco_line %>% filter(Treatment == "Light 25°C"), 
    deco_line %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween ribbons
deco_ribbon_ani <- bind_rows(
  tween_path(
    deco_ribbon %>% filter(Treatment == "Prior"),
    deco_ribbon %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_path(
    deco_ribbon %>% filter(Treatment == "Light 15°C"),
    deco_ribbon %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_path(
    deco_ribbon %>% filter(Treatment == "Dark 15°C"),
    deco_ribbon %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_path(
    deco_ribbon %>% filter(Treatment == "Light 20°C"),
    deco_ribbon %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_path(
    deco_ribbon %>% filter(Treatment == "Light 25°C"),
    deco_ribbon %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Tween densities
deco_dens_ani <- bind_rows(
  tween_polygon(
    deco_dens %>% filter(Treatment == "Prior"),
    deco_dens %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    deco_dens %>% filter(Treatment == "Light 15°C"),
    deco_dens %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_polygon(
    deco_dens %>% filter(Treatment == "Dark 15°C"),
    deco_dens %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_polygon(
    deco_dens %>% filter(Treatment == "Light 20°C"),
    deco_dens %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_polygon(
    deco_dens %>% filter(Treatment == "Light 25°C"),
    deco_dens %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# Shorten prob_dens_ani
prob_dens_ani_short <- bind_rows(
  tween_polygon(
    prob_dens %>% filter(Treatment == "Prior"),
    prob_dens %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 15°C"),
    prob_dens %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Dark 15°C"),
    prob_dens %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 20°C"),
    prob_dens %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_polygon(
    prob_dens %>% filter(Treatment == "Light 25°C"),
    prob_dens %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150)
) %T>%
  print()

# 2.5 Animated figures ####
# Top panel
( ggplot() +
    geom_hline(yintercept = 1) +
    geom_point(
      data = deco_points_ani,
      aes(x, y, colour = colour, alpha = alpha),
      shape = 16, size = 2.35
    ) +
    geom_line(
      data = deco_bars_ani,
      aes(x, y, colour = colour, 
          alpha = alpha, group = Replicate),
    ) +
    geom_line(
      data = deco_line_ani,
      aes(x, y, colour = colour)
    ) +
    geom_polygon(
      data = deco_ribbon_ani,
      aes(x, y, fill = fill, alpha = alpha, group = alpha)
    ) +
    geom_text(
      data = deco_line_ani %>% distinct(Treatment, colour, .frame),
      aes(x = 0.8, y = 0.08, label = Treatment, colour = colour),
      hjust = 0, vjust = 0, size.unit = "pt", size = 12, family = "Futura"
    ) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_y_continuous(labels = scales::label_number(accuracy = c(1, rep(0.1, 4)))) +
    labs(x = "Detrital age (days)", 
         y = "Detrital mass (*m*/*m*<sub><span style='font-size:8.4pt'>0</span></sub>)") +
    coord_cartesian(xlim = c(0, 120), ylim = c(0, 1.6), expand = F, clip = "off") +
    transition_manual(.frame) +
    mytheme +
    theme(
      axis.title.y = element_markdown() # only specify for y axis
    ) ) %>%
  animate(
    fps = 25, # optimal gif frame rate (4 centiseconds per frame) with reasonable file size
    duration = 5 * 3, # number of transitions (5) * time per transition (3 s) to set speed
    detail = 2, # 25*5*3 is 1/2 of all frames (5*150), so the skipped frames can add detail
    width = 20, height = 6,
    units = "cm", res = 300,
    renderer = gifski_renderer(),
    device = "ragg_png" # necessary for italic!
  ) %>%
  anim_save(
    filename = "deco_top.gif", 
    path = here("Figures", "Animation")
  )

# Bottom panel v1
( ggplot() +
    geom_polygon(
      data = deco_dens_ani,
      aes(x = x, y = y, fill = fill)
    ) +
    scale_fill_identity() +
    # Ensure same spacing as top panel
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs(
      x = "Half-life (days)",
      y = expression("Detrital mass ("*italic("m")*"/"*italic("m")[0]*")")
    ) +
    coord_cartesian(xlim = c(0, 120), ylim = c(0, 0.25), expand = F, clip = "off") +
    transition_manual(.frame) +
    mytheme +
    # These will be hidden (transparent) but ensure identical spacing
    theme(
      axis.title.y = element_text(colour = "transparent"),
      axis.text.y = element_text(colour = "transparent"),
      axis.ticks.y = element_line(colour = "transparent"),
      axis.line.y = element_line(colour = "transparent")
    ) ) %>%
  animate(
    fps = 25, duration = 5 * 3, detail = 2,
    width = 20, height = 4, units = "cm", res = 300, 
    renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "deco_bottom_1.gif", 
    path = here("Figures", "Animation")
  )

# Bottom panel v2
( ggplot() +
    geom_polygon(
      data = prob_dens_ani_short,
      aes(x = x, y = y, fill = fill), alpha = 0.8
    ) +
    geom_polygon(
      data = deco_dens_ani,
      aes(x = x, y = y, fill = fill), alpha = 0.8
    ) +
    scale_fill_identity() +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) +
    labs(
      x = "Half-life (days)",
      y = expression("Detrital mass ("*italic("m")*"/"*italic("m")[0]*")")
    ) +
    coord_cartesian(xlim = c(0, 120), ylim = c(0, 0.25), expand = F, clip = "off") +
    transition_manual(.frame) +
    mytheme +
    theme(
      axis.title.y = element_text(colour = "transparent"),
      axis.text.y = element_text(colour = "transparent"),
      axis.ticks.y = element_line(colour = "transparent"),
      axis.line.y = element_line(colour = "transparent")
    ) ) %>%
  animate(
    fps = 25, duration = 5 * 3, detail = 2,
    width = 20, height = 4, units = "cm", res = 300, 
    renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "deco_bottom_2.gif", 
    path = here("Figures", "Animation")
  )

# 3. NMDS ####
# 3.1 Load data ####
ASV_sample <- here("DNA", "RDS", "ASV_sample.rds") %>%
  read_rds() %T>%
  print()

nmds <- ASV_sample %>%
  select(Days, Treatment, MDS1_abund, MDS2_abund) %>%
  rename(Day = Days, MDS1 = MDS1_abund, MDS2 = MDS2_abund) %>%
  mutate(
    Treatment = if_else(Day == 0, "Baseline", Treatment) %>% 
      fct_relevel("Baseline", "Dark 15°C")
  ) %T>%
  print()

nmds_prior_posterior <- here("DNA", "RDS", "nmds_prior_posterior.rds") %>%
  read_rds() %>%
  mutate(
    Treatment = Treatment %>% fct_relevel("Prior", "Baseline", "Dark 15°C")
  ) %T>%
  print()

nmds_time_prior_posterior <- here("DNA", "RDS", "nmds_time_prior_posterior.rds") %>%
  read_rds() %>%
  mutate(
    Treatment = Treatment %>% fct_relevel("Prior", "Baseline", "Dark 15°C")
  ) %T>%
  print()

# 3.2 Aesthetics ####
# Data points
nmds_points <- nmds %>%
  mutate(
    alpha = 0.5, # Baseline opacity
    colour = case_when( # Colours
      Treatment == "Baseline" ~ "#c3b300",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %>%
  rename(x = MDS1, y = MDS2) %>% # Rename to x and y for consistency
  select(Treatment, x, y, alpha, colour) %T>% # Select relevant variables
  print()

nmds_time_points <- nmds %>%
  mutate(
    alpha = 0.5, # Baseline opacity
    colour = case_when( # Colours
      Treatment == "Baseline" ~ "#c3b300",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %>%
  rename(x = MDS1, y = MDS2) %T>%
  print()

# Densities
nmds_dens_x <- nmds_prior_posterior %>%
  filter(k == "MDS1") %>%
  group_by(Treatment) %>% 
  reframe(
    x = c(-3, density(mu, n = 2^10, from = -3, to = 3, bw = 6 * 0.01)$x, 3),
    y = c(0, density(mu, n = 2^10, from = -3, to = 3, bw = 6 * 0.01)$y, 0)
  ) %>%
  group_by(Treatment) %>%
  mutate( y = y * .3 / ( sum(y) * ( x[3] - x[2] ) ) ) %>%
  ungroup() %T>%
  print()

nmds_dens_y <- nmds_prior_posterior %>%
  filter(k == "MDS2") %>%
  group_by(Treatment) %>%
  reframe( # Reversed x and y positions because this density is on the y axis
    x = c(0, density(mu, n = 2^10, from = -2.5, to = 2.5, bw = 5 * 0.01)$y, 0),
    y = c(-2.5, density(mu, n = 2^10, from = -2.5, to = 2.5, bw = 5 * 0.01)$x, 2.5)
  ) %>%
  group_by(Treatment) %>%
  mutate( x = x * .3 / ( sum(x) * ( y[3] - y[2] ) ) ) %>%
  ungroup() %T>%
  print()

nmds_dens_x %<>%
  mutate(
    alpha = 1, # Baseline opacity
    fill = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Baseline" ~ "#c3b300",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %T>%
  print()

nmds_dens_y %<>%
  mutate(
    alpha = 1, # Baseline opacity
    fill = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Baseline" ~ "#c3b300",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    )
  ) %T>%
  print()

# Bivariate densities
require(ggdensity)
require(densityarea) # this is a convenience package building on ggdensity
nmds_ellipse <- nmds_prior_posterior %>%
  pivot_wider(names_from = k, values_from = mu) %>%
  group_by(Treatment) %>%
  reframe(
    density_polygons(
      MDS1, MDS2, probs = c(.5, .8, .9), 
      range_mult = 0, # no padding because these are posteriors
      n = 800 # increase resolution
    )
  ) %T>%
  print()

nmds_ellipse %<>%
  mutate(
    fill = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Baseline" ~ "#c3b300",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    ),
    alpha = case_when(
      prob == .5 ~ .3,
      prob == .8 ~ .4,
      prob == .9 ~ .5
    )
  ) %>%
  select(Treatment, level_id, MDS1, MDS2, fill, alpha) %>%
  rename(x = MDS1, y = MDS2, id = level_id) %T>%
  print()

nmds_time_ellipse <- nmds_time_prior_posterior %>%
  pivot_wider(names_from = k, values_from = mu) %>%
  group_by(Day, Treatment) %>%
  reframe(
    density_polygons(
      MDS1, MDS2, probs = c(.5, .8, .9), 
      range_mult = 0, n = 800
    )
  ) %T>%
  print()

nmds_time_ellipse %<>%
  mutate(
    fill = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Baseline" ~ "#c3b300",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    ),
    alpha = case_when(
      prob == .5 ~ .3,
      prob == .8 ~ .4,
      prob == .9 ~ .5
    ),
    # Convert Day to number
    Day = Day %>% as.character() %>% as.numeric()
  ) %>%
  select(Day, Treatment, level_id, MDS1, MDS2, fill, alpha) %>%
  rename(x = MDS1, y = MDS2, id = level_id) %T>%
  print()

# Text
nmds_text <- nmds_prior_posterior %>%
  distinct(Treatment) %>%
  mutate(
    alpha = 1, # Baseline opacity
    colour = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Baseline" ~ "#c3b300",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    ),
    label = as.character(Treatment),
    x = -2.9, y = 2.2
  ) %T>%
  print()

nmds_time_text <- nmds_time_prior_posterior %>%
  distinct(Day, Treatment) %>%
  arrange(Day, Treatment) %>%
  mutate(
    alpha = 1,
    colour = case_when(
      Treatment == "Prior" ~ "#b5b8ba",
      Treatment == "Baseline" ~ "#c3b300",
      Treatment == "Dark 15°C" ~ "#2e4a5b",
      Treatment == "Light 15°C" ~ "#6a98b4",
      Treatment == "Light 20°C" ~ "#f5a54a",
      Treatment == "Light 25°C" ~ "#d1750c"
    ),
    label = as.character(Treatment),
    Day = Day %>% as.character() %>% as.numeric(),
    x = -2.9, y = 2.2
  ) %T>%
  print()

# 3.3 Static figures ####
# Across days
ggplot() +
  geom_point(
    data = nmds_points,
    aes(x, y, colour = colour, alpha = alpha),
    shape = 16, size = 2.35
  ) +
  geom_polygon(
    data = nmds_ellipse,
    aes(x, y, fill = fill, alpha = alpha, group = id)
  ) +
  geom_polygon(
    data = nmds_dens_x,
    aes(x, y - 2.5, fill = fill, alpha = alpha)
  ) +
  geom_polygon(
    data = nmds_dens_y,
    aes(x - 3, y, fill = fill, alpha = alpha)
  ) +
  geom_text(
    data = nmds_text,
    aes(x, y, label = label, colour = colour, alpha = alpha),
    hjust = 0, vjust = 0, size.unit = "pt", 
    size = 12, family = "Futura" 
  ) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  facet_grid(~ Treatment, space = "free", scales = "free") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-2.5, 2.5), expand = F, clip = "off") +
  mytheme +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linejoin = "mitre"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm")
  )

# By day
ggplot() +
  geom_point(
    data = nmds_time_points,
    aes(x, y, colour = colour, alpha = alpha),
    shape = 16, size = 2.35
  ) +
  geom_polygon(
    data = nmds_time_ellipse,
    aes(x, y, fill = fill, alpha = alpha, group = id)
  ) +
  geom_text(
    data = nmds_time_text,
    aes(x, y, label = label, colour = colour, alpha = alpha),
    hjust = 0, vjust = 0, size.unit = "pt", 
    size = 12, family = "Futura" 
  ) +
  scale_alpha_identity() +
  scale_colour_identity() +
  scale_fill_identity() +
  facet_grid(Day ~ Treatment, space = "free", scales = "free") +
  coord_cartesian(xlim = c(-3, 3), ylim = c(-2.5, 2.5), expand = F, clip = "off") +
  mytheme +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.border = element_rect(linejoin = "mitre"),
    plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm")
  )

# 3.4 Animation ####
# 3.4.1 Across days ####
# Tween points
nmds_points_ani <- bind_rows( 
  tween_state(
    nmds_points %>% filter(Treatment == "Prior"), 
    nmds_points %>% filter(Treatment == "Baseline"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    nmds_points %>% filter(Treatment == "Baseline"), 
    nmds_points %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_state(
    nmds_points %>% filter(Treatment == "Dark 15°C"), 
    nmds_points %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_state(
    nmds_points %>% filter(Treatment == "Light 15°C"), 
    nmds_points %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_state(
    nmds_points %>% filter(Treatment == "Light 20°C"), 
    nmds_points %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_state( 
    # Because prior is sandwiched here, it needs special consideration 
    nmds_points %>% filter(Treatment == "Light 25°C"),
    # Here I am creating a transparent placeholder point rather than 
    # the empty tibble given by setting Treatment == "Prior".
    nmds_points %>% 
      filter(Treatment == "Light 25°C") %>%
      first() %>%
      mutate(alpha = 0),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_state(
    nmds_points %>% filter(Treatment == "Prior"), 
    nmds_points %>% filter(Treatment != "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_state(
    nmds_points %>% filter(Treatment != "Prior"), 
    nmds_points %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150)
) %T>%
  print()

# Tween densities
nmds_dens_x_ani <- bind_rows(
  tween_polygon(
    nmds_dens_x %>% filter(Treatment == "Prior"),
    nmds_dens_x %>% filter(Treatment == "Baseline"),
    ease = "cubic-in-out", nframes = 100,
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    nmds_dens_x %>% filter(Treatment == "Baseline"),
    nmds_dens_x %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100,
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_polygon(
    nmds_dens_x %>% filter(Treatment == "Dark 15°C"),
    nmds_dens_x %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_polygon(
    nmds_dens_x %>% filter(Treatment == "Light 15°C"),
    nmds_dens_x %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100,
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_polygon(
    nmds_dens_x %>% filter(Treatment == "Light 20°C"),
    nmds_dens_x %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100,
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_polygon(
    nmds_dens_x %>% filter(Treatment == "Light 25°C"),
    nmds_dens_x %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100,
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_polygon(
    nmds_dens_x %>% 
      filter(Treatment == "Prior") %>%
      expand_grid(id = 1:5),
    nmds_dens_x %>% 
      filter(Treatment != "Prior") %>%
      droplevels() %>%
      mutate(
        id = as.integer(Treatment),
        alpha = .6
      ),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_polygon(
    nmds_dens_x %>% 
      filter(Treatment != "Prior") %>%
      droplevels() %>%
      mutate(
        id = as.integer(Treatment),
        alpha = .6
      ),
    nmds_dens_x %>% 
      filter(Treatment == "Prior") %>%
      expand_grid(id = 1:5),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150)
) %T>%
  print()

nmds_dens_y_ani <- bind_rows(
  tween_polygon(
    nmds_dens_y %>% filter(Treatment == "Prior"),
    nmds_dens_y %>% filter(Treatment == "Baseline"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    nmds_dens_y %>% filter(Treatment == "Baseline"),
    nmds_dens_y %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_polygon(
    nmds_dens_y %>% filter(Treatment == "Dark 15°C"),
    nmds_dens_y %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_polygon(
    nmds_dens_y %>% filter(Treatment == "Light 15°C"),
    nmds_dens_y %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_polygon(
    nmds_dens_y %>% filter(Treatment == "Light 20°C"),
    nmds_dens_y %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_polygon(
    nmds_dens_y %>% filter(Treatment == "Light 25°C"),
    nmds_dens_y %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_polygon(
    nmds_dens_y %>% 
      filter(Treatment == "Prior") %>%
      expand_grid(id = 1:5),
    nmds_dens_y %>% 
      filter(Treatment != "Prior") %>%
      droplevels() %>%
      mutate(
        id = as.integer(Treatment),
        alpha = .6
      ),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_polygon(
    nmds_dens_y %>% 
      filter(Treatment != "Prior") %>%
      droplevels() %>%
      mutate(
        id = as.integer(Treatment),
        alpha = .6
      ),
    nmds_dens_y %>% 
      filter(Treatment == "Prior") %>%
      expand_grid(id = 1:5),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150)
) %T>%
  print()

# Tween bivariate densities
nmds_ellipse_ani <- bind_rows(
  tween_polygon(
    nmds_ellipse %>% filter(Treatment == "Prior"),
    nmds_ellipse %>% filter(Treatment == "Baseline"),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    nmds_ellipse %>% filter(Treatment == "Baseline"),
    nmds_ellipse %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_polygon(
    nmds_ellipse %>% filter(Treatment == "Dark 15°C"),
    nmds_ellipse %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_polygon(
    nmds_ellipse %>% filter(Treatment == "Light 15°C"),
    nmds_ellipse %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_polygon(
    nmds_ellipse %>% filter(Treatment == "Light 20°C"),
    nmds_ellipse %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_polygon(
    nmds_ellipse %>% filter(Treatment == "Light 25°C"),
    nmds_ellipse %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_polygon(
    nmds_ellipse %>% 
      filter(Treatment == "Prior") %>%
      expand_grid(rep = 1:5) %>%
      mutate(
        id = interaction(id, rep) %>% 
          as.integer(), # Reduce alpha according to overlap
        alpha = 1 - (1 - alpha)^(1/max(rep))
      ) %>%
      select(-rep),
    nmds_ellipse %>% 
      filter(Treatment != "Prior") %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_polygon(
    nmds_ellipse %>% 
      filter(Treatment != "Prior") %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    nmds_ellipse %>% 
      filter(Treatment == "Prior") %>%
      expand_grid(rep = 1:5) %>%
      mutate(
        id = interaction(id, rep) %>% 
          as.integer(),
        alpha = 1 - (1 - alpha)^(1/max(rep))
      ) %>%
      select(-rep),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150)
) %T>%
  print()

# Tween text
nmds_text_ani <- bind_rows( 
  tween_state(
    nmds_text %>% filter(Treatment == "Prior"), 
    nmds_text %>% filter(Treatment == "Baseline"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    nmds_text %>% filter(Treatment == "Baseline"), 
    nmds_text %>% filter(Treatment == "Dark 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1 * 150),
  tween_state(
    nmds_text %>% filter(Treatment == "Dark 15°C"), 
    nmds_text %>% filter(Treatment == "Light 15°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 2 * 150),
  tween_state(
    nmds_text %>% filter(Treatment == "Light 15°C"), 
    nmds_text %>% filter(Treatment == "Light 20°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 3 * 150),
  tween_state(
    nmds_text %>% filter(Treatment == "Light 20°C"), 
    nmds_text %>% filter(Treatment == "Light 25°C"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 4 * 150),
  tween_state(
    nmds_text %>% filter(Treatment == "Light 25°C"), 
    nmds_text %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 5 * 150),
  tween_state(
    nmds_text %>% filter(Treatment == "Prior"), 
    nmds_text %>% 
      filter(Treatment != "Prior") %>%
      mutate(
        x = case_when(
          Treatment == "Baseline" ~ -2.9,
          Treatment == "Dark 15°C" ~ -2.9,
          Treatment == "Light 15°C" ~ -1.41,
          Treatment == "Light 20°C" ~ 0.08,
          Treatment == "Light 25°C" ~ 1.57
        ),
        y = if_else(Treatment == "Baseline", 2.2, 1.9)
      ),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit 
    # specify enter because one label is turning into five
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 6 * 150),
  tween_state(
    nmds_text %>% 
      filter(Treatment != "Prior") %>%
      mutate(
        x = case_when(
          Treatment == "Baseline" ~ -2.9,
          Treatment == "Dark 15°C" ~ -2.9,
          Treatment == "Light 15°C" ~ -1.41,
          Treatment == "Light 20°C" ~ 0.08,
          Treatment == "Light 25°C" ~ 1.57
        ),
        y = if_else(Treatment == "Baseline", 2.2, 1.9)
      ), 
    nmds_text %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
    # specify exit because five labels are turning into one
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 7 * 150)
) %T>%
  print()

# 3.4.2 By day ####
# Tween points
nmds_time_points_ani <- bind_rows( 
  tween_state( # These treatments are only represented at day 0
    nmds_time_points %>% filter(Treatment == "Prior"),
    nmds_time_points %>% filter(Treatment == "Baseline"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    nmds_time_points %>% filter(Treatment == "Baseline"), 
    nmds_time_points %>% filter(Day == 6), # 6 days is the minimum difference
    # linear ease signifies time transition
    ease = "linear", nframes = 60, # 10 frames per day * 6 days = 60 frames
    enter = enter, exit = exit
    # all future enter/exit should also occur across 60 frames
  ) %>%
    mutate(.frame = .frame + 150),
  tween_state(
    nmds_time_points %>% filter(Day == 6), 
    nmds_time_points %>% filter(Day == 13),
    ease = "linear", nframes = 70,
    enter = enter, exit = exit
  ) %>%
    mutate(
      # adjust alpha fade speed to 60 frames
      alpha = if_else(
        .phase %in% c("enter", "exit"),
        deframe(
          tibble(
            seq(0, 0.5, length.out = 70) %>% round(7), # original fade
            c(rep(0, 70-60), seq(0, 0.5, length.out = 60)) # new fade
          ) # round because tweenr sometimes rounds alpha
        )[alpha %>% round(7) %>% as.character()],
        alpha
      ),
      .frame = .frame + 210
    ),
  tween_state(
    nmds_time_points %>% filter(Day == 13), 
    nmds_time_points %>% filter(Day == 27),
    ease = "linear", nframes = 140,
    enter = enter, exit = exit
  ) %>%
    mutate(
      alpha = if_else(
        .phase %in% c("enter", "exit"),
        deframe(
          tibble(
            seq(0, 0.5, length.out = 140) %>% round(7),
            c(rep(0, 140-60), seq(0, 0.5, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        alpha
      ),
      .frame = .frame + 280
    ),
  tween_state(
    nmds_time_points %>% filter(Day == 27), 
    nmds_time_points %>% filter(Day == 37),
    ease = "linear", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    mutate(
      alpha = if_else(
        .phase %in% c("enter", "exit"),
        deframe(
          tibble(
            seq(0, 0.5, length.out = 100) %>% round(7),
            c(rep(0, 100-60), seq(0, 0.5, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        alpha
      ),
      .frame = .frame + 420
    ),
  tween_state( 
    nmds_time_points %>% filter(Day == 37),
    nmds_time_points %>% filter(Day == 82),
    ease = "linear", nframes = 450,
    enter = enter, exit = exit
  ) %>%
    mutate(
      alpha = if_else(
        .phase %in% c("enter", "exit"),
        deframe(
          tibble(
            seq(0, 0.5, length.out = 450) %>% round(7),
            c(rep(0, 450-60), seq(0, 0.5, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        alpha
      ),
      .frame = .frame + 520
    ),
  tween_state(
    nmds_time_points %>% filter(Day == 82), 
    nmds_time_points %>% filter(Day == 119),
    ease = "linear", nframes = 370,
    enter = enter, exit = exit
  ) %>%
    mutate(
      alpha = if_else(
        .phase %in% c("enter", "exit"),
        deframe(
          tibble(
            seq(0, 0.5, length.out = 370) %>% round(7),
            c(rep(0, 370-60), seq(0, 0.5, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        alpha
      ),
      .frame = .frame + 970
    ),
  tween_state(
    nmds_time_points %>% filter(Day == 119), 
    nmds_time_points %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100,
    enter = enter, exit = exit
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1340)
) %T>%
  print()
# 1490 frames in total

nmds_time_points_ani %>% 
  filter(.phase == "transition") %>% 
  ggplot(aes(.frame, Day)) + 
    geom_abline(slope = .1) +
    geom_point() +
    mytheme
# Linear ease maps frames directly onto days

# Tween bivariate densities
nmds_time_ellipse_ani <- bind_rows(
  tween_polygon(
    nmds_time_ellipse %>% filter(Treatment == "Prior"),
    nmds_time_ellipse %>% filter(Treatment == "Baseline"),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50),
  tween_polygon(
    nmds_time_ellipse %>% 
      filter(Treatment == "Baseline") %>%
      expand_grid(rep = 1:5) %>%
      mutate(
        id = interaction(id, rep) %>% 
          as.integer(), # Reduce alpha according to overlap
        alpha = 1 - (1 - alpha)^(1/max(rep))
      ) %>%
      select(-rep),
    nmds_time_ellipse %>% 
      filter(Day == 6) %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ), # linear ease signifies time transition
    ease = "linear", nframes = 60, id = id # 10 frames per day
  ) %>%
    mutate(.frame = .frame + 150),
  tween_polygon(
    nmds_time_ellipse %>% 
      filter(Day == 6) %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    nmds_time_ellipse %>% 
      filter(Day == 13) %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    ease = "linear", nframes = 70, id = id
  ) %>%
    mutate(.frame = .frame + 210),
  tween_polygon(
    nmds_time_ellipse %>% 
      filter(Day == 13) %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    nmds_time_ellipse %>% 
      filter(Day == 27) %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    ease = "linear", nframes = 140, id = id
  ) %>%
    mutate(.frame = .frame + 280),
  tween_polygon(
    nmds_time_ellipse %>% 
      filter(Day == 27) %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    nmds_time_ellipse %>% 
      filter(Day == 37) %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    ease = "linear", nframes = 100, id = id,
    exit = exit # Light 25°C exits here
  ) %>%
    mutate(
      # increase alpha fadeout speed to 60 frames
      alpha = case_when(
        # baseline alpha depends on id
        id %in% c(3, 6, 9, 12) ~ deframe(
          tibble(
            seq(0, .5, length.out = 100) %>% round(7),
            c(rep(0, 100-60), seq(0, .5, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        id %in% c(2, 5, 8, 11) ~ deframe(
          tibble(
            seq(0, .4, length.out = 100) %>% round(7),
            c(rep(0, 100-60), seq(0, .4, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        id %in% c(1, 4, 7, 10) ~ deframe(
          tibble(
            seq(0, .3, length.out = 100) %>% round(7),
            c(rep(0, 100-60), seq(0, .3, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()]
      ),
      .frame = .frame + 420
    ),
  tween_polygon(
    nmds_time_ellipse %>% 
      filter(Day == 37) %>%
      droplevels() %>%
      mutate(
        id = interaction(id, Treatment) %>% 
          as.integer()
      ),
    nmds_time_ellipse %>% 
      filter(Day == 82) %>% # only Light 15°C is left
      mutate(id = id + 3L), # make id consistent
    ease = "linear", nframes = 450, id = id,
    exit = exit
  ) %>%
    mutate(
      alpha = case_when(
        id %in% c(3, 6, 9, 12) ~ deframe(
          tibble(
            seq(0, .5, length.out = 450) %>% round(7),
            c(rep(0, 450-60), seq(0, .5, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        id %in% c(2, 5, 8, 11) ~ deframe(
          tibble(
            seq(0, .4, length.out = 450) %>% round(7),
            c(rep(0, 450-60), seq(0, .4, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        id %in% c(1, 4, 7, 10) ~ deframe(
          tibble(
            seq(0, .3, length.out = 450) %>% round(7),
            c(rep(0, 450-60), seq(0, .3, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()]
      ),
      .frame = .frame + 520
    ),
  tween_polygon(
    nmds_time_ellipse %>% filter(Day == 82),
    nmds_time_ellipse %>% filter(Day == 119),
    ease = "linear", nframes = 370, id = id
  ) %>%
    mutate(.frame = .frame + 970),
  tween_polygon(
    nmds_time_ellipse %>% filter(Day == 119),
    nmds_time_ellipse %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100, id = id
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1340)
) %T>%
  print()
# 1490 frames in total

# Tween text
nmds_time_text_ani <- bind_rows( 
  tween_state(
    nmds_time_text %>% filter(Treatment == "Prior"), 
    nmds_time_text %>% filter(Treatment == "Baseline"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50),
  tween_state(
    nmds_time_text %>% filter(Treatment == "Baseline"), 
    nmds_time_text %>% 
      filter(Day == 6) %>%
      mutate(
        x = case_when(
          Treatment == "Dark 15°C" ~ -2.9,
          Treatment == "Light 15°C" ~ -1.41,
          Treatment == "Light 20°C" ~ 0.08,
          Treatment == "Light 25°C" ~ 1.57
        )
      ),
    ease = "linear", nframes = 60, enter = enter
    # specify enter because one label is turning into four
  ) %>%
    keep_state(nframes = 210) %>% # Keep until day 27
    mutate(.frame = .frame + 150),
  tween_state(
    nmds_time_text %>% 
      filter(Day == 27) %>%
      mutate(
        x = case_when(
          Treatment == "Dark 15°C" ~ -2.9,
          Treatment == "Light 15°C" ~ -1.41,
          Treatment == "Light 20°C" ~ 0.08,
          Treatment == "Light 25°C" ~ 1.57
        )
      ), 
    nmds_time_text %>% 
      filter(Day == 37) %>%
      mutate(
        x = case_when(
          Treatment == "Dark 15°C" ~ -2.9,
          Treatment == "Light 15°C" ~ -1.41,
          Treatment == "Light 20°C" ~ 0.08,
          Treatment == "Light 25°C" ~ 1.57
        )
      ),
    ease = "linear", nframes = 100, exit = exit
    # specify exit because Light 25°C leaves
  ) %>%
    mutate(
      alpha = if_else(
        .phase == "exit",
        deframe(
          tibble( # baseline alpha = 1
            seq(0, 1, length.out = 100) %>% round(7),
            c(rep(0, 100-60), seq(0, 1, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        alpha
      ),
      .frame = .frame + 420
    ),
  tween_state(
    nmds_time_text %>% 
      filter(Day == 37) %>%
      # Re-arrange so rownumbers/ids are paired correctly
      arrange(Treatment %>% fct_relevel("Light 15°C")) %>%
      mutate(
        x = case_when(
          Treatment == "Dark 15°C" ~ -2.9,
          Treatment == "Light 15°C" ~ -1.41,
          Treatment == "Light 20°C" ~ 0.08,
          Treatment == "Light 25°C" ~ 1.57
        )
      ), 
    nmds_time_text %>% 
      filter(Day == 82) %>% # only Light 15°C is left
      mutate(x = -1.41),
    ease = "linear", nframes = 450, exit = exit
    # specify exit because only Light 15°C is left at the end
  ) %>%
    keep_state(nframes = 370) %>% # keep Light 15°C until day 119
    mutate(
      alpha = if_else(
        .phase == "exit",
        deframe(
          tibble(
            seq(0, 1, length.out = 450) %>% round(7),
            c(rep(0, 450-60), seq(0, 1, length.out = 60))
          )
        )[alpha %>% round(7) %>% as.character()],
        alpha
      ),
      .frame = .frame + 520
    ),
  tween_state(
    nmds_time_text %>% 
      filter(Day == 119) %>%
      mutate(x = -1.41), 
    nmds_time_text %>% filter(Treatment == "Prior"),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1340)
) %T>%
  print()
# 1490 frames in total

# Simplify for Day label
nmds_time_day <- nmds_time_text %>% 
  mutate(y = 1.9) %>%
  distinct(Day, x, y)

nmds_time_day_ani <- bind_rows(
  keep_state(
    nmds_time_day %>% filter(Day == 0),
    nframes = 150
  ),
  tween_state(
    nmds_time_day %>% filter(Day == 0), 
    nmds_time_day %>% filter(Day == 119),
    ease = "linear", nframes = 1190
  ) %>%
    mutate(.frame = .frame + 150),
  tween_state(
    nmds_time_day %>% filter(Day == 119), 
    nmds_time_day %>% filter(Day == 0),
    ease = "cubic-in-out", nframes = 100
  ) %>%
    keep_state(nframes = 50) %>%
    mutate(.frame = .frame + 1340)
) %T>%
  print()
# 1490 frames in total

# 2.5 Animated figures ####
# Across days
( ggplot() +
    geom_point(
      data = nmds_points_ani,
      aes(x, y, colour = colour, alpha = alpha),
      shape = 16, size = 2.35
    ) +
    geom_polygon(
      data = nmds_ellipse_ani,
      aes(x, y, fill = fill, alpha = alpha, group = id)
    ) +
    geom_polygon(
      data = nmds_dens_x_ani,
      aes(x, y - 2.5, fill = fill, alpha = alpha)
    ) +
    geom_polygon(
      data = nmds_dens_y_ani,
      aes(x - 3, y, fill = fill, alpha = alpha)
    ) +
    geom_text(
      data = nmds_text_ani,
      aes(x, y, label = label, colour = colour, alpha = alpha),
      hjust = 0, vjust = 0, size.unit = "pt",
      size = 12, family = "Futura"
    ) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    coord_cartesian(
      xlim = c(-3, 3), ylim = c(-2.5, 2.5), 
      expand = F, clip = "off"
    ) +
    transition_manual(.frame) +
    mytheme +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(linejoin = "mitre"),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm")
    ) ) %>%
  animate(
    fps = 25, duration = 8 * 3, detail = 2,
    width = 10, height = 10, units = "cm", res = 300,
    renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "nmds.gif", 
    path = here("Figures", "Animation")
  )

# By day
require(glue)
( ggplot() +
    geom_point(
      data = nmds_time_points_ani,
      aes(x, y, colour = colour, alpha = alpha),
      shape = 16, size = 2.35
    ) +
    geom_polygon(
      data = nmds_time_ellipse_ani,
      aes(x, y, fill = fill, alpha = alpha, group = id)
    ) +
    geom_text(
      data = nmds_time_text_ani,
      aes(x, y, label = label, colour = colour, alpha = alpha),
      hjust = 0, vjust = 0, size.unit = "pt",
      size = 12, family = "Futura"
    ) +
    geom_text(
      data = nmds_time_day_ani, # consistently round half days up 
      aes(x, y, label = glue("Day {floor(Day + .5)}")),
      hjust = 0, vjust = 0, size.unit = "pt",
      size = 12, family = "Futura"
    ) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    coord_cartesian(
      xlim = c(-3, 3), ylim = c(-2.5, 2.5), 
      expand = F, clip = "off"
    ) +
    transition_manual(.frame) +
    mytheme +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(linejoin = "mitre"),
      plot.margin = margin(0.2, 0.2, 0.2, 0.2, unit = "cm")
    ) ) %>%
  animate( # 150 frames per "transition", so 1490/150 "transitions"
    fps = 25, duration = 1490 / 150 * 3, detail = 2,
    width = 10, height = 10, units = "cm", res = 300,
    renderer = gifski_renderer(), device = "ragg_png"
  ) %>%
  anim_save(
    filename = "nmds_time.gif", 
    path = here("Figures", "Animation")
  )
