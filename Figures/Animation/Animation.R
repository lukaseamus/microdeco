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

# Predict across predictor range
source("functions.R")
require(tidybayes)
prob_prediction <- prob_prior_posterior %>%
  spread_continuous(
    data = prob,
    predictor_name = "Day"
  ) %>%
  mutate(
    P_mu = plogis( -5 / mu * ( Day - mu ) ),
    GP_P = rbeta( n() , P_mu * nu , (1 - P_mu) * nu )
  ) %>%
  group_by(Day, Treatment) %>%
  median_qi(P_mu, GP_P, .width = c(.5, .8, .9)) %T>%
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
  mutate( y = y * 2.6 / ( sum(y) * ( x[3] - x[2] ) ) ) %>%
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
    aes(x, y, fill = fill, alpha = alpha, 
        # Interaction grouping is only needed in the static version
        group = interaction(fill, alpha))
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
  facet_grid(~ Treatment, space = "free", scales = "free") +
  labs(x = "Half-life (days)") +
  coord_cartesian(xlim = c(0, 120), expand = F, clip = "off") +
  mytheme +
  theme(axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

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
require(ggtext) # for Markdown syntax in gif

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
    )
  ) %>%
  animate(
    nframes = 750, fps = 25,
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
    geom_point(data = grazing_jitter_ani,
               aes(x = x, y = y, alpha = alpha),
               shape = 16, size = 2.5, colour = "#7030a5") +
    geom_polygon(data = grazing_dens_ani,
                 aes(x = x, y = y, fill = fill)) +
    scale_alpha_identity() +
    scale_fill_identity() +
    coord_cartesian(xlim = c(0, 100), ylim = c(-1, 2), 
                    expand = FALSE, clip = "off") +
    xlab("Defecation (%)<sup><span style='color:white;font-size:8.4pt'>−1</span></sup>") +
    transition_manual(.frame) +
    mytheme +
    theme(axis.title = element_markdown(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank()) ) %>%
  animate(nframes = 450, duration = 15,
          width = 21 * 1/2, height = 10,
          units = "cm", res = 300, renderer = gifski_renderer()) %>%
  anim_save(filename = "grazing_right.gif", path = here("Figures", "Animations"))



# 2. Decomposition ####



# 3. NMDS ####


