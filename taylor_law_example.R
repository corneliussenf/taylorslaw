
library(tidyverse)
library(patchwork)

n <- 100000
m1 <- 0.5
m2 <- 1.5
v0 <- m1 / 2
v1 <- v0 * m2^(0)
v2 <- v0 * m2^(1)
v3 <- v0 * m2^(2)
d0 <- rnorm(n, m1, v0)
d1 <- rnorm(n, m2, v1)
d2 <- rnorm(n, m2, v2)
d3 <- rnorm(n, m2, v3)

p1 <- ggplot() +
  geom_line(stat = "density", aes(x = d0),
            adjust = 3,
            trim = TRUE,
            col = "#000000") +
  geom_line(stat = "density", aes(x = d1),
            adjust = 3,
            trim = TRUE,
            col = "#BBBBBB") +
  xlim(-1.5, 4) +
  geom_vline(xintercept = mean(d0), col = "#000000", linetype = "dotted") +
  geom_vline(xintercept = mean(d1), col = "#BBBBBB", linetype = "dotted") +
  #geom_vline(xintercept = max(d1), col = "#EE6677", linetype = "dashed") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        #axis.line.y = element_blank(),
        axis.ticks = element_blank())
  #annotate("text", x = -1, y = 0.8, label = bquote(b==0), size = 3) +
  #annotate("text", x = max(d1), y = 0.05, label = "*", col = "#EE6677", size = 8) +
  #annotate("text", x = max(d1) + 0.3, y = 0.45, label = "Maximum", col = "#EE6677", size = 2.5)

p2 <- ggplot() +
  geom_line(stat = "density", aes(x = d0), 
            adjust = 3,
            trim = TRUE,
            col = "#000000") +
  geom_line(stat = "density", aes(x = d2), 
            adjust = 3, 
            col = "#BBBBBB",
            trim = TRUE) +
  xlim(-1.5, 4) +
  geom_vline(xintercept = mean(d0), col = "#000000", linetype = "dotted") +
  geom_vline(xintercept = mean(d2), col = "#BBBBBB", linetype = "dotted") +
  #geom_vline(xintercept = max(d2), col = "#EE6677", linetype = "dashed") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        #axis.line.y = element_blank(),
        axis.ticks = element_blank())
  #annotate("text", x = -1, y = 0.8, label = bquote(b==1), size = 3) +
  #annotate("text", x = max(d2), y = 0.05, label = "*", col = "#EE6677", size = 8)

p3 <- ggplot() +
  geom_line(stat = "density", aes(x = d0), 
            adjust = 3,
            trim = TRUE,
            col = "#000000") +
  geom_line(stat = "density", aes(x = d3), 
            adjust = 3,
            trim = TRUE, 
            col = "#BBBBBB") +
  xlim(-1.5, 4) +
  geom_vline(xintercept = mean(d0), col = "#000000", linetype = "dotted") +
  geom_vline(xintercept = mean(d3), col = "#BBBBBB", linetype = "dotted") +
  #geom_vline(xintercept = max(d3), col = "#EE6677", linetype = "dashed") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        #axis.line.y = element_blank(),
        axis.ticks = element_blank()) 
  # annotate("text", x = -1, y = 0.8, label = bquote(b==2), size = 3) +
  # annotate("text", x = max(d3), y = 0.05, label = "*", col = "#EE6677", size = 8)

p1 + p2 + p3 + plot_layout(ncol = 1)

ggsave("Desktop/taylor_law_example.png", width = 3, height = 3)




f = 24
f_c = 3
t = 1 / f;

t = seq(0, 5, t)

A = t

carrier = sin(2*pi*f_c*t)

out = A * carrier

p11 <- ggplot() +
  geom_line(aes(x = t, y = carrier + seq(0, 9, length.out = length(out))),
            col = "#000000") +
  geom_line(aes(x = t, y = seq(0, 9, length.out = length(out))),
            col = "#BBBBBB") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 6)) +
  labs(title = bquote("a) No change in variance")) +
  annotate("text", x = max(t), y = -2.5, label = "time", size = 2) +
  annotate("text", x = -0.45, y = 10.5, label = "disturbance activity", size = 2, angle = 90) +
  coord_cartesian(xlim = c(0, 5), ylim = c(-1, 14), clip = "off")

p22 <- ggplot() +
  geom_line(aes(x = t, y = out + seq(0, 9, length.out = length(out))),
            col = "#000000") +
  geom_line(aes(x = t, y = seq(0, 9, length.out = length(out))),
            col = "#BBBBBB") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 6)) +
  labs(title = "b) Variance changes with mean") +
  annotate("text", x = max(t), y = -2.5, label = "time", size = 2) +
  annotate("text", x = -0.45, y = 10.5, label = "disturbance activity", size = 2, angle = 90) +
  coord_cartesian(xlim = c(0, 5), ylim = c(-1, 14), clip = "off") 

p11 + p22 + 
  (p1 +
  annotate("text", x = 3.1, y = -0.16, label = "disturbance activity", size = 2) +
  coord_cartesian(ylim = c(0, 1.75), clip = "off")) + 
  (p3 +
  annotate("text", x = 3.1, y = -0.16, label = "disturbance activitiy", size = 2) +
  coord_cartesian(ylim = c(0, 1.75), clip = "off")) + 
  plot_layout(ncol = 2)

ggsave("figures/taylor_law_example.png", width = 3.5, height = 3)




