
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
    (params_tmp$sd)^2 * r^2 + 
      2 * (params_tmp$sd) * (params_tmp$mn) * r +
      (params_tmp$mn)^2
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
  
  # ggplot() +
  #   geom_line(
  #     stat = "density",
  #     data = draws,
  #     aes(
  #       x = sqrt(value),
  #       group = key
  #     ),
  #     alpha = 0.3
  #   ) +
  #   geom_density(
  #       data = dat_annual_tmp %>%
  #         filter(
  #           id == i
  #         ),
  #       aes(
  #         x = (rate)
  #       ),
  #       col = "darkred",
  #       linewidth = 1.2
  #     )
  
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

ggplot() +
  geom_violin(
    data = draws,
    aes(
      x = factor(id),
      y = sqrt(value)
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
  )

ggsave("figures/squareroot_assumption.png",
       width = 7.5, 
       height = 4)

