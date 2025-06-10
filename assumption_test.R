
# Packages ----

library(tidyverse)
library(patchwork)
library(sf)
library(ggh4x)

# Colors ----

colors_agents <- c("#3182bd", "#de2d26", "#31a354")
names_agents <- c("Bark beetle & wind", "Wildfire", "Human")

names_subbiomes <- c("Temperte broadleaved/mixed",
                     "Temperate coniferous",
                     "Boreal forest",
                     "Temperate grasslands",
                     "Tundra",
                     "Mediterranean")

# Load data ----

dat_annual <- read_csv("temp/dat_annual_sqrt_test.csv")

ggplot(
  data = dat_annual %>%
    filter(
      id %in% sample(unique(dat_annual$id), 100)
    )
) +
  geom_line(
    stat = "density",
    aes(
      x = sqrt(rate),
      group = id
    ),
    alpha = 0.3,
    adjust = 2
  )

params <- dat_annual %>%
  group_by(
    id,
    agent_grouped
  ) %>%
  summarize(
    mn = mean(rate),
    sd = sd(rate)
  )

# Compare srt- and log-normal distribution ----

### Square-root

n <- 10

draws_out <- vector("list", n)

dat_annual_tmp <- dat_annual %>%
  split(list(.$id, .$agent_grouped)) %>%
  discard(.p = function(x) nrow(x) <= 30) %>%
  bind_rows() %>%
  split(.$agent_grouped) %>%
  map(
    .,
    ~ filter(., id %in% sample(unique(.$id), n))
  ) %>%
  bind_rows() %>%
  split(list(.$id, .$agent_grouped)) %>%
  discard(.p = function(x) nrow(x) == 0)

for (k in 1:length(dat_annual_tmp)) {
  
  d <- dat_annual_tmp[[k]]
  
  params_tmp <- params %>%
    filter(
      id == unique(d$id) & agent_grouped == unique(d$agent_grouped)
    )
  
  fun <- function(x) {
    r <- rnorm(1)
    sd <- params_tmp$sd
    mn <- params_tmp$mn
    rr <- (sd)^2 * r^2 + 2 * (sd) * (mn) * r + (mn)^2
    return(sqrt(rr))
  }
  
  draws <- replicate(
    n = 10000,
    expr = fun()
  )
  
  draws_out[[k]] <- draws %>%
    as_tibble() %>%
    gather() %>%
    mutate(
      value = value
    ) %>%
    mutate(
      id = unique(d$id),
      agent_grouped = unique(d$agent_grouped)
    )
  
}

draws <- draws_out %>%
  bind_rows() %>%
  mutate(
    agent_grouped = case_when(
      agent_grouped == "harvest" ~ "Human",
      agent_grouped == "fire" ~ "Wildfire",
      agent_grouped == "barkbeetle_windthrow" ~ "Bark beetle & wind"
    ),
    agent_grouped = ordered(agent_grouped, levels = c("Bark beetle & wind", "Wildfire", "Human"))
  )

p1 <- ggplot() +
  geom_violin(
    data = draws,
    aes(
      x = factor(id),
      y = (value)
    ),
    adjust = 2
  ) +
  geom_point(
    data = dat_annual_tmp %>%
      bind_rows() %>%
      mutate(
        agent_grouped = case_when(
          agent_grouped == "harvest" ~ "Human",
          agent_grouped == "fire" ~ "Wildfire",
          agent_grouped == "barkbeetle_windthrow" ~ "Bark beetle & wind"
        ),
        agent_grouped = ordered(agent_grouped, levels = c("Bark beetle & wind", "Wildfire", "Human"))
      ),
    aes(
      x = factor(id),
      y = rate,
      col = agent_grouped
    ),
    position = position_jitter(width = 0.05),
    alpha = 0.2
  ) +
  #scale_y_log10() +
  theme_classic() +
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
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 3, name = "Set1")[c(2, 1, 3)]) +
  coord_flip() +
  facet_wrap(~agent_grouped, scales = "free") +
  labs(
    y = bquote("Annual disturbance rate (% "*yr.^-1*")"),
    x = "Randomly selected grid cell"
  ) +
  scale_y_log10()

ggsave("figures/squareroot_assumption.png",
       p1,
       width = 7.5, 
       height = 4)

### Log-normal

n <- 10

draws_out <- vector("list", n)

dat_annual_tmp <- dat_annual %>%
  split(list(.$id, .$agent_grouped)) %>%
  discard(.p = function(x) nrow(x) <= 30) %>%
  bind_rows() %>%
  split(.$agent_grouped) %>%
  map(
    .,
    ~ filter(., id %in% sample(unique(.$id), n))
  ) %>%
  bind_rows() %>%
  split(list(.$id, .$agent_grouped)) %>%
  discard(.p = function(x) nrow(x) == 0)

for (k in 1:length(dat_annual_tmp)) {
  
  d <- dat_annual_tmp[[k]]
  
  params_tmp <- params %>%
    filter(
      id == unique(d$id) & agent_grouped == unique(d$agent_grouped)
    )
  
  fun <- function(x) {
    r <- rnorm(1)
    sd <- params_tmp$sd
    mn <- params_tmp$mn
    mmn <- log(mn) - 0.5 * log((sd/mn)^2 + 1)
    ssd <- sqrt(log((sd/mn)^2 + 1))
    rr <- r * ssd + mmn
    return(exp(rr))
  }
  
  draws <- replicate(
    n = 10000,
    expr = fun()
  )
  
  draws_out[[k]] <- draws %>%
    as_tibble() %>%
    gather() %>%
    mutate(
      value = value
    ) %>%
    mutate(
      id = unique(d$id),
      agent_grouped = unique(d$agent_grouped)
    )
  
}

draws <- draws_out %>%
  bind_rows() %>%
  mutate(
    agent_grouped = case_when(
      agent_grouped == "harvest" ~ "Human",
      agent_grouped == "fire" ~ "Wildfire",
      agent_grouped == "barkbeetle_windthrow" ~ "Bark beetle & wind"
    ),
    agent_grouped = ordered(agent_grouped, levels = c("Bark beetle & wind", "Wildfire", "Human"))
  )

p2 <- ggplot() +
  geom_violin(
    data = draws,
    aes(
      x = factor(id),
      y = (value)
    ),
    adjust = 2
  ) +
  geom_point(
    data = dat_annual_tmp %>%
      bind_rows() %>%
      mutate(
        agent_grouped = case_when(
          agent_grouped == "harvest" ~ "Human",
          agent_grouped == "fire" ~ "Wildfire",
          agent_grouped == "barkbeetle_windthrow" ~ "Bark beetle & wind"
        ),
        agent_grouped = ordered(agent_grouped, levels = c("Bark beetle & wind", "Wildfire", "Human"))
      ),
    aes(
      x = factor(id),
      y = rate,
      col = agent_grouped
    ),
    position = position_jitter(width = 0.05),
    alpha = 0.2
  ) +
  #scale_y_log10() +
  theme_classic() +
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
  scale_color_manual(values = RColorBrewer::brewer.pal(n = 3, name = "Set1")[c(2, 1, 3)]) +
  coord_flip() +
  facet_wrap(~agent_grouped, scales = "free") +
  labs(
    y = bquote("Annual disturbance rate (% "*yr.^-1*")"),
    x = "Randomly selected grid cell"
  ) +
  scale_y_log10()

ggsave("figures/lognormal_assumption.png",
       p2,
       width = 7.5, 
       height = 4)

p1 + p2 + plot_layout(ncol = 1)

ggsave("figures/distribution_assumption.png",
       width = 7.5, 
       height = 8)

# Test for power law ----

resolutions <- list.files("data", pattern = "disturbance_summary*") %>%
  gsub("disturbance_summary_res", "", .) %>%
  gsub(".csv", "", .)

dat <- list.files("data", pattern = "disturbance_summary*", full.names = TRUE) %>%
  map(read_csv) %>%
  set_names(resolutions) %>%
  bind_rows(.id = "resolution") %>%
  mutate(resolution = as.integer(resolution))

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

ggplot(
  data = dat_mod %>%
    mutate(
      x_temp = mean * 100,
      y_temp = var * 100
    )
  ) +
  geom_point(
    aes(
      x = x_temp,
      y = y_temp,
      col = agent
    )
  ) +
  geom_smooth(
    aes(
      x = x_temp,
      y = y_temp,
      col = agent
    ),
    method = "lm"
  )

dat_mod_tmp <- dat_mod %>%
  mutate(
    x_lin = mean * 100,
    x_log = log10(mean * 100),
    y_lin = var * 100,
    y_log = log10(var * 100)
  ) %>%
  select(
    agent,
    resolution,
    x_lin,
    x_log,
    y_lin,
    y_log
  )

newdata <- data.frame(
  x_lin = seq(0.0001, 50, length.out = 100),
  x_log = log10(seq(0.0001, 50, length.out = 100)),
  y_lin = seq(0.0001, 50, length.out = 100),
  y_log = log10(seq(0.0001, 50, length.out = 100))
)

fits_linear <- dat_mod_tmp %>%
  split(list(.$resolution, .$agent)) %>%
  map(., ~ lm(
    y_lin ~ x_lin,
    data = .
  ))
predictions_fits_linear <- fits_linear %>%
  map(
    ~ cbind(
      newdata,
      pred = predict(., newdata = newdata)
    )
  ) %>%
  bind_rows(
    .id = "id"
  ) %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution))
r2_fits_linear <- fits_linear %>%
  map(broom::glance) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  select(agent, resolution, r.squared) %>%
  spread(
    key = resolution,
    value = r.squared
  )

fits_loglinear <- dat_mod_tmp %>%
  split(list(.$resolution, .$agent)) %>%
  map(., ~ lm(
    y_log ~ x_lin,
    data = .
  ))
predictions_fits_loglinear <- fits_loglinear %>%
  map(
    ~ cbind(
      newdata,
      pred = predict(., newdata = newdata)
    )
  ) %>%
  bind_rows(
    .id = "id"
  ) %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution))
r2_fits_loglinear <- fits_loglinear %>%
  map(broom::glance) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  select(agent, resolution, r.squared) %>%
  spread(
    key = resolution,
    value = r.squared
  )

fits_loglog <- dat_mod_tmp %>%
  split(list(.$resolution, .$agent)) %>%
  map(., ~ lm(
    y_log ~ x_log,
    data = .
  ))
predictions_fits_loglog <- fits_loglog %>%
  map(
    ~ cbind(
      newdata,
      pred = predict(., newdata = newdata)
    )
  ) %>%
  bind_rows(
    .id = "id"
  ) %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution))
r2_fits_loglog <- fits_loglog %>%
  map(broom::glance) %>%
  bind_rows(.id = "id") %>%
  separate("id", c("resolution", "agent"), "\\.") %>%
  mutate(resolution = as.integer(resolution)) %>%
  select(agent, resolution, r.squared) %>%
  spread(
    key = resolution,
    value = r.squared
  )

agent.labs <- c("Bark Beetle & Wind", "Wildfire", "Human")
names(agent.labs) <- c("barkbeetle_windthrow", "fire", "harvest")

list(
  r2_fits_linear,
  r2_fits_loglinear,
  r2_fits_loglog
) %>%
  set_names(
    c("Linear", "Log-Linear", "Log-Log")
  ) %>%
  bind_rows(
    .id = "model"
  ) %>%
  pivot_longer(
    cols = `10000`:`160000`,
    names_to = "resolution",
    names_transform = as.integer,
    values_to = "r2"
  ) %>%
  ggplot(
    data = .
  ) +
  geom_bar(
    aes(
      x = model,
      y = r2,
      fill = agent
    ),
    stat = "identity",
    position = "dodge"
  ) +
  facet_grid(
    resolution~agent, 
    labeller = labeller(agent = agent.labs)
  ) +
  theme_classic() +
  scale_fill_manual(values = colors_agents, labels = names_agents) +
  labs(x = "Model",
       y = "Coefficient of determination",
       col = NULL,
       fill = NULL
  ) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    axis.title = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(angle = 90, hjust = 0.5, size = 8)
  ) +
  guides(
    x = "axis_truncated",
    y = "axis_truncated"
  )

ggsave("figures/model_assumption.png",
       width = 5.5, 
       height = 8)


