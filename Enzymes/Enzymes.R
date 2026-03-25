require(here)
require(tidyverse)
require(magrittr)

enzymes <- here("Enzymes", "Enzymes.csv") %>%
  read_csv() %>%
  mutate(
    Aminopeptidase = if_else(Aminopeptidase < 0, 1e-3, Aminopeptidase),
    Glucosidase = if_else(Glucosidase < 0, 1e-3, Glucosidase),
    BG_LAP = Glucosidase / Aminopeptidase, # stoichiometry
    Treatment = case_when(
      Light == 0 ~ "Dark 15°C",
      Temperature == 15 ~ "Light 15°C",
      Temperature == 20 ~ "Light 20°C",
      Temperature == 25 ~ "Light 25°C"
    )
  )%T>%
  print()

# Summarise across technical replicates
enzymes_summary <- enzymes %>%
  summarise(
    across(c(Aminopeptidase, Glucosidase, BG_LAP), list(mean = mean, sd = sd)),
    .by = c(ID, Day, Treatment, Tank)
  ) %T>%
  print()

# Assign baseline samples random treatments
set.seed(1)
random <- enzymes_summary %>%
  filter(is.na(Treatment)) %>%
  distinct(ID) %>%
  mutate(
    Treatment = enzymes_summary %$% Treatment %>% na.omit() %>% 
      unique() %>% rep(3) %>% c(sample(., 3))
    # Temperature = Treatment %>% str_extract("\\d+") %>% as.numeric(),
    # PAR = if_else(Treatment %>% str_detect("Dark"), 0, 8)
  ) %>%
  mutate(
    Tank = case_when(
      Treatment == "Light 15°C" ~ sample(1:4, n()),
      Treatment == "Dark 15°C" ~ sample(5:8, n()),
      Treatment == "Light 20°C" ~ sample(9:12, n()),
      Treatment == "Light 25°C" ~ sample(13:16, n())
    ),
    .by = Treatment
  ) %T>%
  print()

# Add randomly generated treatments to enzymes_summary
enzymes_summary %<>%
  left_join(random, by = "ID") %>%
  mutate(
    Treatment = coalesce(Treatment.x, Treatment.y),
    Tank = coalesce(Tank.x, Tank.y)
    # Temperature = coalesce(Temperature.x, Temperature.y),
    # PAR = coalesce(PAR.x, PAR.y)
  ) %>%
  select(-c(ends_with(".x"), ends_with(".y"))) %T>%
  print()

# Calculate means across replicates for geom_line (same origin)
enzymes_medians <- enzymes_summary %>%
  mutate(Treatment = if_else(Day == 0, "Baseline", Treatment)) %>%
  summarise(
    Glucosidase = median(Glucosidase_mean),
    Aminopeptidase = median(Aminopeptidase_mean),
    BG_LAP = median(BG_LAP_mean),
    .by = c(Day, Treatment)
  ) %>%
  mutate(count = if_else(Day == 0, 4, 1)) %>% uncount(count) %>%
  mutate(
    Treatment = case_when(
      row_number() == 1 ~ "Light 15°C",
      row_number() == 2 ~ "Dark 15°C",
      row_number() == 3 ~ "Light 20°C",
      row_number() == 4 ~ "Light 25°C",
      TRUE ~ Treatment
    )
  ) %T>%
  print()


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

require(ggh4x)
Fig_2a <- ggplot() +
    geom_point(
      data = enzymes_summary,
      aes(Day, Glucosidase_mean, colour = Treatment),
      size = 2.5, shape = 16, alpha = 0.7
    ) +
    # geom_pointrange(
    #   data = enzymes_summary,
    #   aes(Day, Glucosidase_mean,
    #       ymin = Glucosidase_mean - Glucosidase_sd,
    #       ymax = Glucosidase_mean + Glucosidase_sd,
    #       colour = Treatment),
    #   shape = 16, alpha = 0.5
    # ) +
    geom_line(
      data = enzymes_medians, 
      aes(Day, Glucosidase, colour = Treatment),
      lineend = "round", linejoin = "round"
    ) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                        guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                      guide = "none") +
    scale_y_continuous(breaks = seq(0, 6, 2)) +
    facet_grid(~Treatment, scales = "free", space = "free") +
    facetted_pos_scales(
      x = list(
        Treatment == "Light 15°C" ~
          scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)),
        Treatment != "Light 15°C" ~
          scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 30))
      )
    ) +
    labs(x = "Detrital age (days)",
         y = expression("BG (µg cm"^-2*" h"^-1*")")) +
    coord_cartesian(ylim = c(0, 6), expand = F, clip = "off") +
    mytheme +
    theme(axis.title.y = element_text(vjust = 0, margin = margin(t = 0)))
    
Fig_2a

Fig_2b <- ggplot() +
    geom_point(
      data = enzymes_summary,
      aes(Day, Aminopeptidase_mean, colour = Treatment),
      size = 2.5, shape = 16, alpha = 0.7
    ) +
    geom_line(
      data = enzymes_medians, 
      aes(Day, Aminopeptidase, colour = Treatment),
      lineend = "round", linejoin = "round"
    ) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                        guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                      guide = "none") +
    facet_grid(~Treatment, scales = "free", space = "free") +
    facetted_pos_scales(
      x = list(
        Treatment == "Light 15°C" ~
          scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)),
        Treatment != "Light 15°C" ~
          scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 30))
      )
    ) +
    labs(x = "Detrital age (days)",
         y = expression("LAP (µg cm"^-2*" h"^-1*")")) +
    coord_cartesian(ylim = c(0, 80), expand = F, clip = "off") + 
    mytheme +
    theme(axis.title.y = element_text(vjust = 0, margin = margin(t = 0)))
    
Fig_2b

Fig_2c <- ggplot() +
    geom_point(
      data = enzymes_summary,
      aes(Day, BG_LAP_mean, colour = Treatment),
      size = 2.5, shape = 16, alpha = 0.7
    ) +
    geom_line(
      data = enzymes_medians, 
      aes(Day, BG_LAP, colour = Treatment),
      lineend = "round", linejoin = "round"
    ) +
    geom_hline(yintercept = 1) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                        guide = "none") +
    scale_fill_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                      guide = "none") +
    scale_y_continuous(labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))) +
    facet_grid(~Treatment, scales = "free", space = "free") +
    facetted_pos_scales(
      x = list(
        Treatment == "Light 15°C" ~
          scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)),
        Treatment != "Light 15°C" ~
          scale_x_continuous(limits = c(0, 60), breaks = seq(0, 60, 30))
      )
    ) +
    labs(x = "Detrital age (days)",
         y = "BG / LAP") +
    coord_cartesian(ylim = c(0, 1.5), expand = F, clip = "off") + 
    mytheme +
    theme(axis.title.y = element_text(vjust = 0.16))
    
Fig_2c

require(patchwork)
Fig_2 <- ( 
  ( Fig_2a +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.margin = margin(0, 0.5, 0.2, 0, unit = "cm")) ) / 
  ( Fig_2b +
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            strip.text = element_blank(),
            plot.margin = margin(0.5, 0.5, 0.2, 0, unit = "cm")) ) / 
  ( Fig_2c +
      theme(strip.text = element_blank(),
            plot.margin = margin(0.5, 0.5, 0.2, 0, unit = "cm")) )
)
  # plot_annotation(tag_levels = c("a", "b", "c")) &
  # theme(plot.tag = element_text(family = "Futura", size = 15, face = "bold"))
Fig_2

Fig_2 %>%
  ggsave(filename = "Fig_2.pdf", path = "Figures",
         device = cairo_pdf, height = 15, width = 20, units = "cm")



# Load carbon and nitrogen data
deco_summary <- here("Decomposition", "RDS", "deco_summary.rds") %>%
  read_rds %T>%
  print()

# Join data
C_N <- deco_summary %>%
  left_join(
    enzymes_summary %>%
      mutate(Treatment = Treatment %>% fct(),
             Tank = Tank %>% factor())
  ) %>%
  rename(BG_LAP = BG_LAP_mean) %>%
  drop_na(C_N, BG_LAP) %T>%
  print()

Fig_3 <- C_N %>%
  left_join(
    C_N %>% summarise(
      BG_LAP_median = median(BG_LAP),
      C_N_median = median(C_N),
      .by = c(Day, Treatment)
    )
  ) %>%
  ggplot() +
    geom_point(
      aes(C_N, BG_LAP, colour = Treatment),
      size = 2.5, shape = 16, alpha = 0.2
    ) +
    geom_segment(
      aes(x = C_N, xend = C_N_median,
          y = BG_LAP, yend = BG_LAP_median,
          colour = Treatment),
      alpha = 0.2, lineend = "round"
    ) +
    geom_path(
      data = . %>% distinct(Day, Treatment, C_N_median, BG_LAP_median),
      aes(C_N_median, BG_LAP_median, colour = Treatment),
      arrow = arrow(length = unit(0.2, "cm"),
                    type = "closed", angle = 20),
      lineend = "round", linejoin = "round"
    ) +
    geom_hline(yintercept = 1) +
    # geom_vline(xintercept = 1) +
    scale_colour_manual(values = c("#2e4a5b", "#6a98b4", "#f5a54a", "#d1750c"), 
                        guide = "none") +
    #scale_y_continuous(labels = scales::label_number(accuracy = c(1, 0.1, 1, 0.1))) +
    scale_x_log10(breaks = c(10, 20, 40, 80)) +
    scale_y_log10(breaks = c(0.01, 0.1, 1, 10),
                  labels = scales::label_number(accuracy = c(0.01, 0.1, 1, 1))) +
    facet_grid(~Treatment) +
    labs(x = "C / N", y = "BG / LAP") +
    coord_cartesian(ylim = c(0.01, 10), xlim = c(10, 80),
                    expand = F, clip = "off") +
    mytheme

Fig_3

Fig_3 %>%
  ggsave(filename = "Fig_3.pdf", path = "Figures",
         device = cairo_pdf, height = 6, width = 20, units = "cm")
  