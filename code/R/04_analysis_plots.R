# Load necessary libraries
library("ggplot2")
library("cowplot")
ggplot2::theme_set(theme_cowplot())
library("gridExtra")
library("dplyr")
library("tidyr")
library("here")
library("optparse")

source(here::here("R", "00_utils.R"))

# define arguments -------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, "--bnab", default = "10E8+PG9+VRC07",
                     help = "The bnAb combination of interest")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
print(args)

bnabs <- unlist(strsplit(args$bnab, "+", fixed = TRUE))
plot_bnab <- paste0(bnabs, collapse = "_")
nice_bnab_combo <- janitor::make_clean_names(args$bnab)
nice_bnabs <- janitor::make_clean_names(bnabs)

plots_dir <- here::here("..", "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}

# analysis dataset
analysis_dataset <- readRDS(file = here::here("..", "data", paste0("analysis_dataset_", nice_bnab_combo, ".rds")))

# plots of neutralization sensitivity values (either measured or combined) -----

# number of variables:
ncol(analysis_dataset %>% 
       select(-starts_with("seq"),
              -starts_with("ic"),
              -starts_with("sens_")))
# pivot to longer format: one row per virus/outcome/antibody
long_analysis_dataset <- analysis_dataset %>% 
  select(seq.id.lanl, starts_with("ic"), starts_with("sens_")) %>% 
  pivot_longer(cols = -seq.id.lanl,
               names_to = c("outcome", "bnab"),
               names_pattern = "(ic50|ic80|sens_ic50|sens_ic80)_(.*)",
               values_to = "neutralization") %>% 
  mutate(combo_method = case_when(
    bnab == "additive" ~ "additive",
    bnab == "bh" ~ "Bliss-Hill",
    (bnab != "additive") & (bnab != "bh") ~ "none"
  ), bnab = case_when(
    combo_method == "additive" ~ nice_bnab_combo,
    combo_method == "Bliss-Hill" ~ nice_bnab_combo,
    combo_method == "none" ~ bnab
  ), bnab = factor(plot_nab_name(bnab, bnabs = c(bnabs, args$bnab),
                                 nice_bnabs = c(nice_bnabs, plot_bnab)), 
                   levels = c(bnabs, args$bnab)))

continuous_value_boxplot <- long_analysis_dataset %>% 
  filter(outcome != "sens_ic50", outcome != "sens_ic80") %>% 
  ggplot(aes(x = bnab, y = neutralization, color = combo_method)) +
  geom_boxplot() +
  labs(y = "Neutralization value", x = "Broadly neutralizing antibody regimen",
       color = "Combination method") +
  facet_grid(cols = vars(outcome)) +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
  

continuous_value_histogram <- long_analysis_dataset %>% 
  filter(outcome != "sens_ic50", outcome != "sens_ic80") %>% 
  ggplot(aes(x = neutralization)) +
  geom_histogram() +
  labs(x = "Neutralization value") +
  facet_grid(cols = vars(bnab), rows = vars(outcome, combo_method),
             labeller = labeller(bnab = label_both, outcome = label_parsed, combo_method = label_parsed))

ggsave(filename = paste0(plots_dir, "/ic50ic80_boxplot_", nice_bnab_combo, ".png"),
       continuous_value_boxplot, width = 8, height = 4, units = "in", dpi = 300)

ggsave(filename = paste0(plots_dir, "/ic50ic80_histogram_", nice_bnab_combo, ".png"),
       continuous_value_histogram, width = 8, height = 8, units = "in", dpi = 300)


# numbers susceptible etc. for IC50, IC80 --------------------------------------

binary_outcome_table <- long_analysis_dataset %>%
  group_by(bnab, combo_method) %>%
  summarize(
    ic50_percent = round(sum(neutralization[outcome == "sens_ic50" & !is.na(neutralization)]) / sum(!is.na(neutralization[outcome == "sens_ic50"])) * 100, 2),
    ic50_size = sum(!is.na(neutralization[outcome == "sens_ic50"])),
    ic80_percent = round(sum(neutralization[outcome == "sens_ic80" & !is.na(neutralization)]) / sum(!is.na(neutralization[outcome == "sens_ic80"])) * 100, 2),
    ic80_size = sum(!is.na(neutralization[outcome == "sens_ic80"]))
  ) %>%
  mutate(across(ends_with("_size"), ~ ifelse(. == 0, NA, .))) %>%
  na_if(NA) %>%
  rename(
    `% susceptible (IC50)` = ic50_percent,
    `number of viruses (IC50)` = ic50_size,
    `% susceptible (IC80)` = ic80_percent,
    `number of viruses (IC80)` = ic80_size
  )
write.csv(binary_outcome_table, file = paste0(plots_dir, "/binary_outcome_table_", nice_bnab_combo, ".csv"))

# plot prediction performance --------------------------------------------------
outcome_order <- c("true_combo_ic50", "ic50", "true_combo_ic80", "ic80", 
                   "true_combo_sens_ic50", "sens_ic50", "true_combo_sens_ic80", "sens_ic80")
nice_outcomes <- c("Observed IC50", "IC50", "Observed IC80", "IC80", 
                   "Observed susceptibility (IC50)", "Susceptibility (IC50)", 
                   "Observed susceptibility (IC80)", "Susceptibility (IC80)")
pred_perf <- readRDS(file = here::here("..", "results", paste0("perf_", nice_bnab_combo, ".rds"))) |> 
  mutate(bnab = factor(plot_nab_name(bnab, bnabs = c(bnabs, args$bnab),
                                     nice_bnabs = c(nice_bnabs, plot_bnab)), 
                       levels = c(bnabs, args$bnab)),
         combo_time = factor(combo_time, levels = c("before", "after"),
                             labels = c("Combine then predict (CP)", "Predict then combine (PC)")),
         outcome = factor(outcome, levels = outcome_order, labels = nice_outcomes, ordered = TRUE)) |> 
  filter(!(grepl("true", outcome) & grepl("individual", method)))

# Mapping for colors and shapes
color_ba_mapping <- c("Combine then predict (CP)" = "lightblue", "Predict then combine (PC)" = "coral3")
shape_method_mapping <- c("individual" = "circle", "additive" = "triangle", "Bliss-Hill" = "square")

###### Continuous ######
########################
continuous_bnab_perf <- pred_perf %>%
  filter(!grepl("Susc", outcome, ignore.case = TRUE))

# Create the plot with switched "before" and "after" positions and custom x-axis order
cont_ind_plot <- continuous_bnab_perf |> 
  filter(method == "individual") |> 
  ggplot(aes(x = outcome, y = perf, shape = bnab)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), position = position_dodge(width = 0.2), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = "Outcome", y = "CV R-squared", shape = "bnAb") +
  ggtitle("Individual bnAbs") +
  theme(legend.position = "bottom")
cont_comb_plot <- continuous_bnab_perf |> 
  filter(method != "individual") |> 
  ggplot(aes(x = outcome, y = perf, color = combo_time, shape = method)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), position = position_dodge(width = 0.2), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  scale_shape_manual(values = shape_method_mapping) +
  scale_color_manual(values = color_ba_mapping) +
  labs(x = "Outcome", y = "CV R-squared", color = "Approach", shape = "Method") +
  ggtitle("bnAb combination") +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))

full_cont_plot <- cowplot::plot_grid(cont_ind_plot, cont_comb_plot, ncol = 2)
# Show the plot
ggsave(filename = paste0(plots_dir, "/continuous_performance_", nice_bnab_combo, ".png"), 
       full_cont_plot, width = 10, height = 6, dpi = 300)

###### Binary ######
########################

binary_bnab_perf <- pred_perf %>%
  filter(grepl("Susc", outcome, ignore.case = TRUE))

# Create the plot with switched "before" and "after" positions and custom x-axis order
bin_ind_plot <- binary_bnab_perf |> 
  filter(method == "individual") |> 
  ggplot(aes(x = outcome, y = perf, shape = bnab)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), position = position_dodge(width = 0.2), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  labs(x = "Outcome", y = "CV AUC", shape = "bnAb") +
  ggtitle("Individual bnAbs") +
  theme(legend.position = "bottom")
bin_comb_plot <- binary_bnab_perf |> 
  filter(method != "individual") |> 
  ggplot(aes(x = outcome, y = perf, color = combo_time, shape = method)) +
  geom_point(position = position_dodge(width = 0.2), size = 3) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), position = position_dodge(width = 0.2), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
  scale_shape_manual(values = shape_method_mapping) +
  scale_color_manual(values = color_ba_mapping) +
  labs(x = "Outcome", y = "CV AUC", color = "Approach", shape = "Method") +
  ggtitle("bnAb combination") +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom") +
  scale_x_discrete(guide = guide_axis(n.dodge = 2))
full_bin_plot <- cowplot::plot_grid(bin_ind_plot, bin_comb_plot, ncol = 2)
# Show the plot
ggsave(filename = paste0(plots_dir, "/binary_performance_", nice_bnab_combo, ".png"), 
       full_bin_plot, width = 10, height = 6, dpi = 300)

# combined performance plot
ind_legend <- cowplot::get_legend(cont_ind_plot)
comb_legend <- cowplot::get_legend(cont_comb_plot)
all_perf_plot <- cowplot::plot_grid(
  cowplot::plot_grid(cont_ind_plot + theme(legend.position = "none"), 
                     bin_ind_plot + theme(legend.position = "none") + ggtitle(""), 
                     ind_legend, nrow = 3, rel_heights = c(1, 1, .15)),
  cowplot::plot_grid(cont_comb_plot + theme(legend.position = "none"), 
                     bin_comb_plot + theme(legend.position = "none") + ggtitle(""), 
                     comb_legend, nrow = 3, rel_heights = c(1, 1, .15)),
  ncol = 2, rel_widths = c(0.7, 1)
)
ggsave(filename = paste0(plots_dir, "/all_performance_", nice_bnab_combo, ".png"),
       all_perf_plot, width = 11.5, height = 6, dpi = 300)
