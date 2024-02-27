
# Packages ----

library(tidyverse)
library(patchwork)
library(sf)
library(ggh4x)

# Colors ----

#colors_agents <- wesanderson::wes_palette("Zissou1", n = 5)[c(1, 5, 3)]
colors_agents <- c("#3182bd", "#de2d26", "#31a354")
names_agents <- c("Bark beetle & wind", "Wildfire", "Human")

names_subbiomes <- c("Temperte broadleaved/mixed",
                     "Temperate coniferous",
                     "Boreal forest",
                     "Temperate grasslands",
                     "Tundra",
                     "Mediterranean")

# colors_biomes <- wesanderson::wes_palette("Chevalier1")[c(1, 2, 3)]
# names_biomes <- c("Boreal",
#                   "Mediterranean",
#                   "Temperate")

# Load data ----

resolutions <- list.files("data", pattern = "disturbance_summary*") %>%
  gsub("disturbance_summary_res", "", .) %>%
  gsub(".csv", "", .)

dat <- list.files("data", pattern = "disturbance_summary*", full.names = TRUE) %>%
  map(read_csv) %>%
  set_names(resolutions) %>%
  bind_rows(.id = "resolution") %>%
  mutate(resolution = as.integer(resolution))

head(dat)

summary(dat)

ecoregions <- list.files("data", pattern = "ecoregion_res*", full.names = TRUE) %>%
  map(read_csv) %>%
  set_names(resolutions) %>%
  bind_rows(.id = "resolution") %>%
  mutate(resolution = as.integer(resolution))

dat_mod <- dat %>%
  rename(agent = agent_grouped) %>%
  mutate(
    log_mean_sqrt = log10(mean_sqrt * 100),
    log_var_sqrt = log10(var_sqrt * 100),
    log_mean = log10(mean * 100),
    log_var = log10(var * 100)
  ) %>%
  na.omit() %>%
  filter(var > 0) %>%
  left_join(ecoregions, by = c("resolution", "id")) %>%
  mutate(
    subbiome = factor(biome, labels = names_subbiomes),
    biome = case_when(
      biome %in% c(4, 5, 8) ~ "Temperate",
      biome %in% c(6, 11) ~ "Boreal",
      biome == 12 ~ "Mediterranean"
    )
  )

# Check biomes ----

res <- 40000

grid <- read_sf(paste0("data/grid_", res, ".gpkg"))

grid_tmp <- grid %>%
  full_join(
    x = .,
    y = dat_mod %>%
      filter(resolution == res) %>%
      group_by(id) %>%
      summarize(
        biome = unique(biome, na.rm = TRUE)
      )
  )

ggplot() +
  geom_sf(
    data = grid_tmp %>% filter(!is.na(biome)),
    aes(
      fill = biome,
      col = biome
    )
  ) +
  theme_void() +
  labs(
    fill = "Biome",
    col = "Biome"
  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    )
  )

ggsave(
  filename = "figures/map_biome.png", 
  width = 5, 
  height = 5
  )

# Plot data ----
  
grid_plot <- grid %>%
  right_join(
    dat_mod %>%
      filter(resolution == res),
    multiple = "all",
    by = c("id")
  ) %>%
  mutate(
    mean = ifelse(is.na(mean), 0, mean)
  ) %>%
  mutate(
    agent = case_when(
      agent == "harvest" ~ "Human",
      agent == "fire" ~ "Wildfire",
      agent == "barkbeetle_windthrow" ~ "Bark beetle & wind"
    ),
    agent = ordered(agent, levels = c("Bark beetle & wind", "Wildfire", "Human"))
  )

p_mean <- ggplot() +
    geom_sf(
      data = grid %>%
        filter(id %in% unique(grid_plot$id)),
      fill = "grey",
      col = NA
    ) +
    geom_sf(
    data = grid_plot,
    aes(
      fill = log10(mean)
    ),
    col = NA
  ) +
  facet_wrap(
    ~ agent
  ) +
  theme_void() +
  labs(
    fill = bquote("log"[10]*"%"),
    col = NULL
  ) +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.2, "cm"),
    plot.title = element_text(size = 10)
  ) +
  scale_fill_viridis_c()

p_var <- ggplot() +
  geom_sf(
    data = grid %>%
      filter(id %in% unique(grid_plot$id)),
    fill = "grey",
    col = NA
  ) +
  geom_sf(
    data = grid_plot,
    aes(
      fill = log10(var)
    ),
    col = NA
  ) +
  facet_wrap(
    ~ agent
  ) +
  theme_void() +
  labs(
    fill = bquote("log"[10]*"%"),
    col = NULL
  ) +
  theme(
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.2, "cm"),
    plot.title = element_text(size = 10)
  ) +
  scale_fill_viridis_c()

p_mean +
  p_var +
  plot_layout(ncol = 1) +
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  )

ggsave(
  "figures/mean_variance_maps.png",
  width = 7.5,
  height = 5.5
  )

# Estimate power-law coefficients ----

fits <- dat_mod %>%
  split(list(.$resolution, .$agent)) %>%
  map(., ~ lmodel2::lmodel2(
    log_var ~ log_mean, 
    data = .,
    nperm = 99,
    range.x = "interval",
    range.y = "interval"
  ))

fits_biome <- dat_mod %>%
  split(list(.$resolution, .$agent, .$biome)) %>%
  discard(.p = function(x) nrow(x) < 5) %>%
  map(., ~ lmodel2::lmodel2(
    log_var ~ log_mean, 
    data = .,
    nperm = 99,
    range.x = "interval",
    range.y = "interval"
  ))

# Model performance ----

r2 <- fits %>%
  map(broom::glance) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  select(agent, resolution, r.squared) %>%
  spread(
    key = resolution,
    value = r.squared
  )
  
write_csv(
  x = r2,
  file = "figures/model_r2.csv"
  )

r2_biome <- fits_biome %>%
  map(broom::glance) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent", "biome"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  select(biome, agent, resolution, r.squared) %>%
  spread(
    key = resolution,
    value = r.squared
  )

write_csv(
  r2_biome,
  file = "figures/model_r2_biome.csv"
)

fits %>%
  map(broom::glance) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  ggplot(data = .) +
  geom_point(aes(x = resolution, 
                 y = r.squared, 
                 col = agent)) +
  geom_line(aes(x = resolution, 
                y = r.squared, 
                col = agent, 
                group = agent)) +
  theme_classic() +
  scale_color_manual(values = colors_agents, labels = names_agents) +
  scale_fill_manual(values = colors_agents, labels = names_agents) +
  labs(x = parse(text = "'Spatial grain ('*10^4*' m)'"),
       y = "Coefficient of determination",
       col = NULL,
       fill = NULL
  ) +
  scale_x_log10(
    breaks = sort(as.integer(resolutions)),
    labels = c(1, 2, 4, 8, 16)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    limits = c(0, 1)
  ) +
  theme(
    legend.position = c(0, 0),
    legend.justification = c(-0.05, -0.1),
    legend.background = element_blank(),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.1, "cm"),
    legend.text = element_text(size = 7),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8)
  ) +
  guides(
    x = "axis_truncated",
    y = "axis_truncated"
  )

ggsave(
  filename = "figures/model_r2_spatial_scaling.png", 
  width = 3, 
  height = 3
  )

fits_biome %>%
  map(broom::glance) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent", "biome"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  mutate(
    agent = factor(agent, labels = names_agents)
  ) %>%
  ggplot(data = .) +
  geom_point(aes(x = resolution, 
                 y = r.squared, 
                 col = agent)) +
  geom_line(aes(x = resolution, 
                y = r.squared, 
                col = agent, 
                group = agent)) +
  theme_classic() +
  scale_color_manual(values = colors_agents, labels = names_agents) +
  scale_fill_manual(values = colors_agents, labels = names_agents) +
  labs(x = parse(text = "'Spatial grain ('*10^4*' m)'"),
       y = "Coefficient of determination",
       col = NULL,
       fill = NULL
  ) +
  scale_x_log10(
    breaks = sort(as.integer(resolutions)),
    labels = c(1, 2, 4, 8, 16)
  ) +
  scale_y_continuous(
    breaks = c(0, 0.5, 1),
    limits = c(0, 1)
  ) +
  theme(
    legend.position = c(0, 0),
    legend.justification = c(-0.05, -0.1),
    legend.background = element_blank(),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.1, "cm"),
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8)
  ) +
  guides(
    x = "axis_truncated",
    y = "axis_truncated"
  ) +
  facet_wrap(~biome)

ggsave(
  filename = "figures/model_r2_spatial_scaling_biomes.png", 
  width = 7.5, 
  height = 2.5
)

r2 %>%
  gather(
    key = resolution,
    value = r2,
    -agent
  ) %>%
  group_by(
    agent
  ) %>%
  summarize(
    r2 = mean(r2)
  )

# Variation in coefficients over spatial scales ----

fits %>%
  map(broom::tidy) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  filter(term == "Slope" & method == "RMA") %>%
  group_by(agent) %>%
  summarise(
    m = mean(estimate),
    sd = sd(estimate),
    mn = min(estimate),
    mx = max(estimate)
  )

fits_biome %>%
  map(broom::tidy) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent", "biome"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  filter(term == "Slope" & method == "RMA") %>%
  group_by(agent, biome) %>%
  summarise(
    n = n(),
    m = mean(estimate),
    sd = sd(estimate),
    mn = min(estimate),
    mx = max(estimate)
  ) %>%
  mutate(
    biome = ordered(biome, levels = c("Boreal", "Temperate", "Mediterranean"))
  )

plotdat_all <- fits %>%
  map(broom::tidy) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  filter(term == "Slope" & method == "RMA")
  
ggplot(data = plotdat_all) +
  geom_point(aes(x = resolution, 
                 y = estimate, 
                 col = agent)) +
  geom_line(aes(x = resolution, 
                y = estimate, 
                col = agent, 
                group = agent)) +
  geom_ribbon(aes(x = resolution,
                    ymin = conf.low, 
                    ymax = conf.high,
                    fill = agent,
                    group = agent), 
              alpha = 0.3,
              show.legend = FALSE
              ) +
  theme_classic() +
  scale_color_manual(values = colors_agents, labels = names_agents) +
  scale_fill_manual(values = colors_agents, labels = names_agents) +
  labs(x = parse(text = "'Spatial grain ('*10^4*' m)'"),
       y = "Power law coefficient",
       col = NULL,
       fill = NULL
       ) +
  scale_x_log10(
    breaks = sort(as.integer(resolutions)),
    labels = c(1, 2, 4, 8, 16)
  ) +
  scale_y_continuous(
    breaks = c(0, 1, 2, 3, 4, 5),
    limits = c(0, 5)
  ) +
  theme(
    legend.position = c(0, 0),
    legend.justification = c(-0.05, -0.1),
    legend.background = element_blank(),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.1, "cm"),
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8)
  ) +
  guides(
    x = "axis_truncated",
    y = "axis_truncated"
    )

ggsave(
  filename = "figures/power_law_coefficient_spatial_scaling.png", 
  width = 2.5, 
  height = 2.5
  )

plotdat_biome <- fits_biome %>%
  map(broom::tidy) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent", "biome"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  filter(term == "Slope" & method == "RMA")

plotdat_biome %>%
  mutate(biome = ordered(biome, levels = c("Boreal", "Temperate", "Mediterranean"))) %>%
  ggplot(.) +
  geom_point(aes(x = resolution, 
                 y = estimate, 
                 col = agent)) +
  geom_line(aes(x = resolution, 
                y = estimate, 
                col = agent, 
                group = agent)) +
  geom_ribbon(aes(x = resolution,
                  ymin = conf.low, 
                  ymax = conf.high,
                  fill = agent,
                  group = agent), 
              alpha = 0.3,
              show.legend = FALSE
  ) +
  theme_classic() +
  scale_color_manual(values = colors_agents, labels = names_agents) +
  scale_fill_manual(values = colors_agents, labels = names_agents) +
  labs(
    x = parse(text = "'Spatial grain ('*10^4*' m)'"),
    #x = "Spatial grain (m)",
    y = "Power law coefficient",
    col = NULL,
    fill = NULL
  ) +
  scale_x_log10(
    breaks = sort(as.integer(resolutions)),
    labels = c(1, 2, 4, 8, 16)
  ) +
  scale_y_continuous(
    breaks = c(-2, -1, 0, 1, 2, 3, 4, 5, 6),
    limits = c(-2, 6)
  ) +
  theme(
    legend.position = c(0, 0),
    legend.justification = c(-0.05, -0.1),
    legend.background = element_blank(),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.1, "cm"),
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8)
  ) +
  # geom_hline(yintercept = c(1, 2),
  #            linetype = "dashed",
  #            col = "#BBBBBB") +
  guides(
    x = "axis_truncated",
    y = "axis_truncated"
  ) +
  facet_wrap(
    ~ biome, 
    scales = "free"
    )

ggsave(
  filename = "figures/power_law_coefficient_spatial_scaling_biome.png", 
  width = 7.5, 
  height = 3
  )

# Power law relationship ----

coefs_bestscale <- fits %>%
  map(broom::tidy) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  filter(
    method == "RMA"
  ) %>%
  left_join(
    x = r2 %>%
      pivot_longer(
        cols = `10000`:`160000`,
        names_to = "resolution",
        names_transform = as.integer,
        values_to = "r2"
      )
  ) %>%
  group_by(
    agent
  ) %>%
  filter(
    r2 == max(r2)
  ) %>%
  ungroup()

write_csv(
  x = coefs_bestscale, 
  file = "figures/coefficients_bestscale.csv"
)

coefs_bestscale_biome <- fits_biome %>%
  map(broom::tidy) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent", "biome"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  filter(
    method == "RMA"
  ) %>%
  left_join(
    x = r2_biome %>%
      pivot_longer(
        cols = `10000`:`160000`,
        names_to = "resolution",
        names_transform = as.integer,
        values_to = "r2"
      )
  ) %>%
  group_by(
    agent,
    biome
  ) %>%
  filter(
    r2 == max(r2)
  ) %>%
  ungroup()

write_csv(
  x = coefs_bestscale_biome, 
  file = "figures/coefficients_bestscale_biome.csv"
)

coefs_bestscale_plot <- fits %>%
  map(broom::tidy) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  filter(method == "RMA") %>%
  filter(
    (agent %in% c("harvest") & resolution == 10000) |
      (agent %in% c("barkbeetle_windthrow", "fire") & resolution == 160000)
  ) %>%
  select(agent, method, estimate, term) %>%
  pivot_wider(
    names_from = "term",
    values_from = "estimate"
  ) %>%
  mutate(
    agent = case_when(
      agent == "harvest" ~ "Human",
      agent == "fire" ~ "Wildfire",
      agent == "barkbeetle_windthrow" ~ "Bark beetle & wind"
    ),
    agent = ordered(agent, levels = c("Bark beetle & wind", "Wildfire", "Human"))
  )

pred <- dat_mod %>%
  filter(
    (agent %in% c("harvest") & resolution == 10000) |
      (agent %in% c("barkbeetle_windthrow", "fire") & resolution == 160000)
  ) %>%
  group_by(agent) %>%
  summarise(
    min_x = min(log_mean),
    max_x = max(log_mean)
  ) %>%
  split(.$agent) %>%
  map2(
    .x = .,
    .y = coefs_bestscale_plot %>%
      split(.$agent),
    ~ data.frame(
      x = seq(.x$min_x, .x$max_x, length.out = 100),
      y = .y$Intercept + .y$Slope * seq(.x$min_x, .x$max_x, length.out = 100)
      )
  ) %>%
  bind_rows(.id = "agent")%>%
  mutate(
    agent = case_when(
      agent == "harvest" ~ "Human",
      agent == "fire" ~ "Wildfire",
      agent == "barkbeetle_windthrow" ~ "Bark beetle & wind"
    ),
    agent = ordered(agent, levels = c("Bark beetle & wind", "Wildfire", "Human"))
  )

ggplot(dat_mod %>%
         filter(
           (agent %in% c("harvest") & resolution == 10000) |
             (agent %in% c("barkbeetle_windthrow", "fire") & resolution == 160000)
         ) %>%
         mutate(
           agent = case_when(
             agent == "harvest" ~ "Human",
             agent == "fire" ~ "Wildfire",
             agent == "barkbeetle_windthrow" ~ "Bark beetle & wind"
           ),
           agent = ordered(agent, levels = c("Bark beetle & wind", "Wildfire", "Human"))
         ),
       aes(
         x = mean * 100,
         y = var * 100,
         col = agent,
         fill = agent
       )
  ) +
  geom_point(
    alpha = 0.4,
    shape = 16
  ) +
  geom_line(
    data = pred,
    aes(
      x = 10^x,
      y = 10^y
    ),
    col = "grey20",
    linewidth = 1
  ) +
  theme_classic() +
  scale_color_manual(values = colors_agents,
                     labels = names_agents) +
  scale_fill_manual(values = colors_agents,
                    labels = names_agents) +
  labs(
    x = bquote("Mean disturbance rate (% "*yr.^-1*")"),
    y = bquote("Temporal variance"),
    col = NULL,
    fill = NULL
  ) +
  theme(
    legend.position = "none",
    legend.justification = c(-0.05, 1.1),
    legend.background = element_blank(),
    legend.key.height = unit(0.1, "cm"),
    legend.key.width = unit(0.1, "cm"),
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8),
    plot.margin = unit(c(0.5, 0, 0, 0.5), "cm")
  ) +
  guides(
    x = "axis_truncated",
    y = "axis_truncated"
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  facet_wrap(
    ~ agent,
    scales = "free_y"
    ) +
  geom_text(
    data = r2 %>%
      pivot_longer(
        cols = `10000`:`160000`
      ) %>%
      filter(
        (agent %in% c("harvest") & name == 10000) |
          (agent %in% c("barkbeetle_windthrow", "fire") & name == 160000)
      ) %>%
      mutate(
        agent = case_when(
          agent == "harvest" ~ "Human",
          agent == "fire" ~ "Wildfire",
          agent == "barkbeetle_windthrow" ~ "Bark beetle & wind"
        ),
        agent = ordered(agent, levels = c("Bark beetle & wind", "Wildfire", "Human"))
      ),
    aes(
      x = 10^-2.5,
      y = 10^1.3,
      label = paste("R^2 ==", round(value, 2))
    ),
    parse = TRUE,
    size = 2.7,
    col = "black"
  ) +
  geom_text(
    data = coefs_bestscale_plot,
    aes(
      x = 10^-2.5,
      y = 10^0,
      label = paste("y==", round(exp(Intercept), 2), "*", "x^", round(Slope, 2))
    ),
    parse = TRUE,
    size = 2.7,
    col = "black"
  ) +
  scale_x_log10(
    breaks = c(10^-4, 10^-2, 10^-0, 10^2),
    limits = c(10^-4.1, 10^2),
    labels = scales::label_log()
  ) +
  scale_y_log10(
    breaks = c(10^-10, 10^-8, 10^-6, 10^-4, 10^-2, 10^0, 10^2),
    limits = c(10^-10, 10^2),
    labels = scales::label_log()
  )
  
ggsave(
  filename = "figures/scatterplot_bestscale_rma.png",
  width = 7,
  height = 2.5
)

# Differences between biomes ----

coefs_biome_bestscale <- fits_biome %>%
  map(broom::tidy) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent", "biome"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  filter(
    term == "Slope" & 
      method == "RMA"
  ) %>%
  left_join(
    x = r2_biome %>%
      pivot_longer(
        cols = `10000`:`160000`,
        names_to = "resolution",
        names_transform = as.integer,
        values_to = "r2"
      )
  ) %>%
  group_by(
    biome,
    agent
  ) %>%
  filter(
    r2 == max(r2)
  ) %>%
  ungroup()

write_csv(
  x = coefs_biome_bestscale, 
  file = "figures/coefficients_bestscale_biome.csv"
)

# Simulations ----

# draws <- function(i, x, a, b) {
#   r <- rnorm(n = i,
#              mean = sqrt(x),
#              sd = sqrt(sqrt(exp(a + b * log(x)))))^2
#   return(r)
# }

draws <- function(i, x, a, b) {
  r <- rnorm(i)
  mn <- x
  sd <- sqrt((exp(a + b * log(x))))
  rr <- (sd)^2 * r^2 + 2 * (sd) * (mn) * r + (mn)^2
  return(rr)
  
  # See here:
  # https://math.stackexchange.com/questions/620045/mean-and-variance-of-squared-gaussian-y-x2-where-x-sim-mathcaln0-sigma
  
}

# f <- function(x, a, b) (exp(a + b * log((x))))
# 
# s <- seq(0.0001, 10, 0.5)
# 
# plot(s, f(s, 0, 0), type = "l", ylim = c(0, 150))
# lines(s, f(s, 0, 1))
# lines(s, f(s, 0, 2))

coefs_simulation <- coefs_bestscale %>%
  filter(agent == "barkbeetle_windthrow")

n <- 100000
a <- coefs_simulation$estimate[1]

sim_res <- expand_grid(
    x = seq(0.1, 3, 0.1),
    b = c(2.3)
  ) %>%
  mutate(id = paste(x, b, sep = "-")) %>%
  split(.$id) %>%
  map(
    ~ data.frame(
        dist_rate = rep(.x$x, n),
        coef = rep(.x$b, n),
        draws = draws(n, .x$x, a, .x$b),
        n = 1:n
      )
  ) %>%
  bind_rows()

simulations_1 <- sim_res %>% 
  filter(dist_rate %in% c(0.5, 1, 2)) %>%
  ggplot(data = .) +
  geom_histogram(
    aes(
      x = draws,
      y = ..density..,
      fill = factor(dist_rate, labels = c("0.5", "1.0", "2.0"))
    )
  ) +
  theme_classic() +
  labs(
    x = bquote("Realized disturbance rates (%"~yr.^-1*")"),
    y = "Density",
    fill = bquote("Mean disturbance rate (%"~yr.^-1*"):")
  ) +
  theme(
    legend.position = c(1, 1),
    legend.justification = c(1, 1.1),
    legend.direction = "horizontal",
    legend.background = element_blank(),
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  ) +
  guides(
    x = "axis_truncated",
    y = "axis_truncated"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1),
      title.position = "top"
      )
    ) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(name = "YlGnBu", n = 7)[c(2, 5, 7)]) +
  xlim(0, 20) +
  ylim(0, 0.6) +
  facet_wrap(
    ~ dist_rate,
    ncol = 1,
    scales = "free_x"
  )

ggsave(
  filename = "figures/simulations.png", 
  plot = simulations_1,
  width = 2.5, 
  height = 2.5
  )

simulations_2 <- sim_res %>%
  group_by(dist_rate, coef) %>%
  summarise(
    pr_2.5 = mean(draws > 2.5),
    pr_5 = mean(draws > 5),
    pr_10 = mean(draws > 10)
  ) %>%
  pivot_longer(
    cols = pr_2.5:pr_10
  ) %>%
  separate(
    col = "name",
    into = c("measure", "threshold"),
    sep = "_"
  ) %>%
  pivot_wider(
    names_from = measure,
    values_from = value
  ) %>%
  ggplot(
    data = .
  ) +
  geom_smooth(
    aes(
      x = dist_rate,
      y = pr,
      col = factor(threshold, 
                   levels = c("2.5", "5", "10"),
                   labels = c("2.5", "5.0", "10.0")
                   )
    ),
    method = "gam",
    se = FALSE
  ) +
  scale_color_manual(values = RColorBrewer::brewer.pal(name = "RdPu", n = 7)[c(2, 5, 7)]) +
  theme_classic() +
  labs(
    x = "Mean disturbance rate (%)",
    y = "Probability",
    col = bquote("Annual disturbance rate (%"~yr.^-1*") greater:")
  ) +
  theme(
    legend.position = c(0, 1),
    legend.justification = c(-0.05, 1.2),
    legend.background = element_blank(),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7)
  ) +
  guides(
    x = "axis_truncated",
    y = "axis_truncated"
  ) +
  xlim(0, 2) +
  ylim(0, 0.8)

ggsave(
  filename = "figures/simulations_probability.png", 
  plot = simulations_2,
  width = 2.5, 
  height = 2.5
  )

simulations_1 + 
  simulations_2 + 
  plot_annotation(
    tag_levels = "a",
    tag_suffix = ")"
  ) +
  plot_layout()

ggsave(
  filename = "figures/simulations_figure.png", 
  width = 7, 
  height = 3.5
)

sim_res %>%
  group_by(dist_rate, coef) %>%
  summarise(
    pr_2.5 = mean(draws > 2.5),
    pr_5 = mean(draws > 5),
    pr_10 = mean(draws > 10)
  ) %>%
  View()
