# create plot for schematic figure

set.seed(1234)
concentrations <- rep(10^seq(-0.75, 2.25, by = 0.25), each = 2)
inhibition <- 1 - pmin(100, pmax(0, plogis(1 - concentrations, location = 0, scale = 0.6) + 
  rnorm(length(concentrations), 0, 0.025)))
example_dat <- tibble::tibble(
  inhibition = round(inhibition * 100, 2),
  concentration = concentrations
)
library("ggplot2")
library("dplyr")
library("cowplot")
theme_set(theme_cowplot())
initial_plot <- example_dat %>% 
  ggplot(aes(x = concentration, y = inhibition)) +
  geom_point(position = position_dodge(width = 0.05), size = 2) +
  scale_x_log10() +
  ylim(c(0, 100)) +
  labs(y = "Percent inhibition", x = "Concentration") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25))
initial_plot

ggsave(here::here("..", "plots", "schematic_neutralization_curve_1.png"),
       initial_plot + 
         geom_smooth(se = FALSE, linewidth = 1.5, span = 0.5) +
         geom_point(position = position_dodge(width = 0.05), size = 2),
       width = 5, height = 5, units = "in", dpi = 300
       )

set.seed(20231030)
concentrations <- rep(10^seq(-0.75, 2.25, by = 0.25), each = 2)
inhibition <- 1 - pmin(100, pmax(0, plogis(1 - concentrations, location = 0, scale = 0.6) + 
                                   rnorm(length(concentrations), 0.05, 0.075)))
example_dat <- tibble::tibble(
  inhibition = round(inhibition * 100, 2),
  concentration = concentrations
)
plot_2 <- example_dat %>% 
  ggplot(aes(x = concentration, y = inhibition)) +
  geom_point(position = position_dodge(width = 0.05), size = 2) +
  scale_x_log10() +
  ylim(c(0, 100)) +
  labs(y = "Percent inhibition", x = "Concentration") +
  theme(axis.text = element_text(size = 20), axis.title = element_text(size = 25))

ggsave(here::here("..", "plots", "schematic_neutralization_curve_2.png"),
       plot_2 + 
         geom_smooth(se = FALSE, linewidth = 1.5, span = 0.7) +
         geom_point(position = position_dodge(width = 0.05), size = 2),
       width = 5, height = 5, units = "in", dpi = 300
)
