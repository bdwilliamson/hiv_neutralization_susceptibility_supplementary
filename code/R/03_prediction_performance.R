# obtain prediction performance for a given bnAb

# load required functions and packages -----------------------------------------
library("here")
library("optparse")
library("cvAUC")
library("vimp")
library("dplyr")

source(here::here("R", "00_utils.R"))

# define arguments -------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, "--bnab", default = "10E8+PG9+VRC07",
                     help = "The bnAb combination of interest")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
print(args)

bnabs <- unlist(strsplit(args$bnab, "+", fixed = TRUE))
nice_bnab_combo <- janitor::make_clean_names(args$bnab)
nice_bnabs <- janitor::make_clean_names(bnabs)

# load the predictions ---------------------------------------------------------
all_preds <- readRDS(file = here::here("..", "results", paste0("preds_", nice_bnab_combo, ".rds"))) |> 
  filter(!(grepl("true_combo_", outcome) & grepl("individual", method)))

# compute prediction performance -----------------------------------------------
continuous_preds <- all_preds %>% 
  filter(!grepl("sens", outcome))

binary_preds <- all_preds %>% 
  filter(grepl("sens", outcome))

continuous_pred_perf <- continuous_preds %>% 
  group_by(bnab, outcome, combo_time, method) %>% 
  summarize(perf = cvR2(predictions = preds, y = y, folds = folds,
                        na.rm = TRUE)$cvR2,
            ci = cvR2(predictions = preds, y = y, folds = folds,
                      na.rm = TRUE)$ci, .groups = "drop")

binary_pred_perf <- binary_preds %>% 
  group_by(bnab, outcome, combo_time, method) %>% 
  summarize(perf = tryCatch(cvAUC(predictions = preds[!is.na(preds)], labels = y[!is.na(preds)], folds = folds[!is.na(preds)])$cvAUC, 
                            error = function(e) cvAUC(predictions = preds[!is.na(preds)], labels = y[!is.na(preds)])$cvAUC),
            ci = tryCatch(matrix(
              ci.cvAUC(predictions = preds[!is.na(preds)], labels = y[!is.na(preds)], folds = folds[!is.na(preds)])$ci,
              nrow = 1
            ), error = function(e) {
              matrix(
                ci.cvAUC(predictions = preds[!is.na(preds)], labels = y[!is.na(preds)])$ci,
                nrow = 1
              )
            }),
            .groups = "drop")

all_pred_perf <- rbind(continuous_pred_perf, binary_pred_perf) %>% 
  mutate(cil = ci[, 1], ciu = ci[, 2]) %>% 
  select(-ci)

# save
saveRDS(all_pred_perf, file = here::here("..", "results", paste0("perf_", nice_bnab_combo, ".rds")))
print("Performance complete.")
