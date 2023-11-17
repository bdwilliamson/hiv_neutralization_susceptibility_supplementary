# get plots, tables of results for combo timing simulation

# load required functions and packages -----------------------------------------
library("dplyr")
library("tidyr")
library("rprojroot")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())

this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)
plots_dir <- paste0(proj_root, "/plots/")

source(paste0(proj_root, "/code/sims/gen_data.R"))
source(paste0(proj_root, "/code/sims/00_utils.R"))

# read in results files --------------------------------------------------------
output_dir <- paste0(proj_root, "/results/sims/")
all_output_files <- list.files(output_dir, pattern = "output")
all_output <- do.call(rbind, lapply(as.list(all_output_files), function(file) readRDS(paste0(output_dir, file))))

# get summaries of prediction performance --------------------------------------
summary_tib <- all_output |> 
  group_by(n, p, model, correlated, strength, outcome, measure, mab) |> 
  mutate(width_init = ciu - cil) |> 
  summarize(mn_point_est = mean(point_est), ESE = sd(point_est), width = mean(width_init),
            .groups = "drop") |> 
  mutate(mab_fct = factor(case_when(
    mab == "ind_1" ~ "Individual bnAb 1",
    mab == "ind_2" ~ "Individual bnAb 2",
    mab == "ind_3" ~ "Individual bnAb 3",
    mab == "comb_pre" ~ "Combine then predict (CP)",
    mab == "comb_post" ~ "Predict then combine (PC)"
  )))

ind_perf_plot_cont <- summary_tib |> 
  filter(outcome == "IC80", grepl("ind", mab, ignore.case = TRUE)) |> 
  ggplot(aes(x = n, y = mn_point_est, shape = mab_fct,
             ymin = mn_point_est - ESE * qnorm(0.975),
             ymax = mn_point_est + ESE * qnorm(0.975))) + 
  geom_errorbar(width = 50, linewidth = 0.5, position = position_dodge(width = 25),
                alpha = 0.5) +
  geom_point(size = 2, position = position_dodge(width = 25), fill = "red") +
  scale_x_continuous(breaks = unique(summary_tib$n),
                   labels = as.character(unique(summary_tib$n)),
                   guide = guide_axis(n.dodge = 2)) +
  scale_shape_manual(values = c(15, 17, 18)) +
  ggtitle("Prediction performance") +
  ylab("CV R-squared") +
  xlab("n") +
  labs(shape = "bnAb") +
  facet_grid(rows = vars(strength))
ind_width_plot_cont <- summary_tib |> 
  filter(outcome == "IC80", grepl("ind", mab, ignore.case = TRUE)) |> 
  ggplot(aes(x = n, y = width, shape = mab_fct)) + 
  geom_point(size = 2, position = position_dodge(width = 25), fill = "red") +
  scale_x_continuous(breaks = unique(summary_tib$n),
                     labels = as.character(unique(summary_tib$n)),
                     guide = guide_axis(n.dodge = 2)) +
  scale_shape_manual(values = c(15, 17, 18)) +
  ggtitle("CI width") +
  ylab("95% CI width") +
  xlab("n") +
  labs(shape = "bnAb") +
  facet_grid(rows = vars(strength))

ind_perf_plot_bin <- summary_tib |> 
  filter(grepl("Susc", outcome, ignore.case = TRUE), grepl("ind", mab, ignore.case = TRUE)) |> 
  ggplot(aes(x = n, y = mn_point_est, shape = mab_fct,
             ymin = mn_point_est - ESE * qnorm(0.975),
             ymax = mn_point_est + ESE * qnorm(0.975))) + 
  geom_errorbar(width = 50, linewidth = 0.5, position = position_dodge(width = 25),
                alpha = 0.5) +
  geom_point(size = 2, position = position_dodge(width = 25), fill = "red") +
  scale_x_continuous(breaks = unique(summary_tib$n),
                     labels = as.character(unique(summary_tib$n)),
                     guide = guide_axis(n.dodge = 2)) +
  scale_shape_manual(values = c(15, 17, 18)) +
  ylab("CV AUC") +
  xlab("n") +
  labs(shape = "bnAb") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  facet_grid(rows = vars(strength))
ind_width_plot_bin <- summary_tib |> 
  filter(outcome != "IC80", grepl("ind", mab, ignore.case = TRUE)) |> 
  ggplot(aes(x = n, y = width, shape = mab_fct)) + 
  geom_point(size = 2, position = position_dodge(width = 25), fill = "red") +
  scale_x_continuous(breaks = unique(summary_tib$n),
                     labels = as.character(unique(summary_tib$n)),
                     guide = guide_axis(n.dodge = 2)) +
  scale_shape_manual(values = c(15, 17, 18)) +
  ylab("95% CI width") +
  xlab("n") +
  labs(shape = "bnAb") +
  facet_grid(rows = vars(strength))

ind_legend <- get_legend(ind_perf_plot_bin)
ind_plot <- cowplot::plot_grid(cowplot::plot_grid(
  plot_grid(ind_perf_plot_cont + theme(legend.position = "none"),
            ind_width_plot_cont + theme(legend.position = "none"),
            ncol = 2),
  plot_grid(ind_perf_plot_bin + theme(legend.position = "none"),
            ind_width_plot_bin + theme(legend.position = "none"),
            ncol = 2),
  nrow = 2
), ind_legend, nrow = 2, rel_heights = c(1, .1))

comb_perf_plot_cont <- summary_tib |> 
  filter(outcome == "IC80", !grepl("ind", mab, ignore.case = TRUE)) |> 
  ggplot(aes(x = n, y = mn_point_est, shape = mab_fct,
             ymin = mn_point_est - ESE * qnorm(0.975),
             ymax = mn_point_est + ESE * qnorm(0.975))) + 
  geom_errorbar(width = 50, linewidth = 0.5, position = position_dodge(width = 25),
                alpha = 0.5) +
  geom_point(size = 2, position = position_dodge(width = 25)) +
  scale_x_continuous(breaks = unique(summary_tib$n),
                     labels = as.character(unique(summary_tib$n)),
                     guide = guide_axis(n.dodge = 2)) +
  scale_shape_manual(values = c(4, 16)) +
  ggtitle("Prediction performance") +
  ylab("CV R-squared") +
  xlab("n") +
  labs(shape = "Approach") +
  facet_grid(rows = vars(strength))
comb_width_plot_cont <- summary_tib |> 
  filter(outcome == "IC80", !grepl("ind", mab, ignore.case = TRUE)) |> 
  ggplot(aes(x = n, y = width, shape = mab_fct)) + 
  geom_point(size = 2, position = position_dodge(width = 25)) +
  scale_x_continuous(breaks = unique(summary_tib$n),
                     labels = as.character(unique(summary_tib$n)),
                     guide = guide_axis(n.dodge = 2)) +
  scale_shape_manual(values = c(4, 16)) +
  ggtitle("CI width") +
  ylab("95% CI width") +
  xlab("n") +
  labs(shape = "Approach") +
  facet_grid(rows = vars(strength))

comb_perf_plot_bin <- summary_tib |> 
  filter(outcome != "IC80", !grepl("ind", mab, ignore.case = TRUE)) |> 
  ggplot(aes(x = n, y = mn_point_est, shape = mab_fct,
             ymin = mn_point_est - ESE * qnorm(0.975),
             ymax = mn_point_est + ESE * qnorm(0.975))) + 
  geom_errorbar(width = 50, linewidth = 0.5, position = position_dodge(width = 25),
                alpha = 0.5) +
  geom_point(size = 2, position = position_dodge(width = 25)) +
  scale_x_continuous(breaks = unique(summary_tib$n),
                     labels = as.character(unique(summary_tib$n)),
                     guide = guide_axis(n.dodge = 2)) +
  scale_shape_manual(values = c(4, 16)) +
  ylab("CV AUC") +
  xlab("n") +
  labs(shape = "Approach") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  facet_grid(rows = vars(strength))
comb_width_plot_bin <-  summary_tib |> 
  filter(outcome != "IC80", !grepl("ind", mab, ignore.case = TRUE)) |> 
  ggplot(aes(x = n, y = width, shape = mab_fct)) + 
  geom_point(size = 2, position = position_dodge(width = 25)) +
  scale_x_continuous(breaks = unique(summary_tib$n),
                     labels = as.character(unique(summary_tib$n)),
                     guide = guide_axis(n.dodge = 2)) +
  scale_shape_manual(values = c(4, 16)) +
  ylab("95% CI width") +
  xlab("n") +
  labs(shape = "Approach") +
  facet_grid(rows = vars(strength))
comb_legend <- get_legend(comb_perf_plot_bin)
comb_plot <- cowplot::plot_grid(cowplot::plot_grid(
  plot_grid(comb_perf_plot_cont + theme(legend.position = "none"),
            comb_width_plot_cont + theme(legend.position = "none"),
            ncol = 2),
  plot_grid(comb_perf_plot_bin + theme(legend.position = "none"),
            comb_width_plot_bin + theme(legend.position = "none"),
            ncol = 2),
  nrow = 2
), comb_legend, nrow = 2, rel_heights = c(1, .1))

ggsave(filename = paste0(plots_dir, "sim_results_individual.png"),
       ind_plot, width = 8, height = 8, units = "in", dpi = 300)
ggsave(filename = paste0(plots_dir, "sim_results_combo.png"),
       comb_plot, width = 8, height = 8, units = "in", dpi = 300)

# get outcome proportion by generating a large dataset
set.seed(20231017)
test_dat <- gen_data(n = 1e6, p = 1000, model = "simple", correlated = FALSE, n_mabs = 3)
all_ic80 <- 10 ^ do.call(cbind, test_dat$y_cont)
ic80_additive <- log10(apply(all_ic80, 1, function(x) additive_combo_method(x)))
susc_additive <- as.numeric(ic80_additive < 0)
cat("Individual susceptibility proportions:\n")
print(apply(test_dat$y_bin, 2, mean))
cat("Combination susceptibility proportion:\n")
print(mean(susc_additive))

set.seed(20231019)
test_dat <- gen_data(n = 1e6, p = 1000, model = "simple", correlated = FALSE, n_mabs = 3,
                     strength = "weak")
all_ic80 <- 10 ^ do.call(cbind, test_dat$y_cont)
ic80_additive <- log10(apply(all_ic80, 1, function(x) additive_combo_method(x)))
susc_additive <- as.numeric(ic80_additive < 0)
cat("Individual susceptibility proportions:\n")
print(apply(test_dat$y_bin, 2, mean))
cat("Combination susceptibility proportion:\n")
print(mean(susc_additive))
