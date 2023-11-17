# get predictions for neutralization susceptibility to a given bnAb

# load required functions and packages -----------------------------------------
library("here")
library("readr")
library("optparse")
library("glmnet")
library("janitor")
library("dplyr")
library("tidyr")
library("parallel")
library("foreach")
library("doParallel")

source(here::here("R", "00_utils.R"))

# define arguments -------------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, "--bnab", default = "10E8+PG9+VRC07",
                     help = "The bnAb combination of interest")
parser <- add_option(parser, "--K", default = 20, 
                     help = "Number of cross-fitting folds")
parser <- add_option(parser, "--parallel", default = TRUE, 
                     help = "Should we run cv.glmnet in parallel?")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
print(args)

print(paste0("Running bnAb ", args$bnab))

# set up the dataset -----------------------------------------------------------
bnabs <- unlist(strsplit(args$bnab, "+", fixed = TRUE))
nice_bnab_combo <- janitor::make_clean_names(args$bnab)
nice_bnabs <- janitor::make_clean_names(bnabs)
bnab_file <- list.files(path = here::here("..", "data"), 
                        pattern = paste0("analysis_dataset_", nice_bnab_combo, "*.rds"))
analysis_dataset <- readRDS(here::here("..", "data", bnab_file))
has_combo_neut <- any(grepl(nice_bnab_combo, names(analysis_dataset)))

n_susceptible <- analysis_dataset %>% 
  select(starts_with("sens"), -contains(nice_bnab_combo)) %>% 
  summarize(across(everything(), .fns = list(is_susceptible = ~ sum(.x, na.rm = TRUE),
                                             non_susceptible = ~ sum(1 - .x, na.rm = TRUE))))
susceptible_vec <- as.numeric(n_susceptible %>% select(contains("is_")))
nonsusceptible_vec <- as.numeric(n_susceptible %>% select(contains("non_")))
binary_K <- min(5, susceptible_vec, nonsusceptible_vec)

# set up specific args
args_ic50 <- c(args, list("outcome" = "ic50"))
args_ic80 <- c(args, list("outcome" = "ic80"))
args_sens_ic50 <- c(args, list("outcome" = "sens_ic50"))
args_sens_ic80 <- c(args, list("outcome" = "sens_ic80"))
args_sens_ic50$K <- binary_K
args_sens_ic80$K <- binary_K
args_ic50_before <- c(args_ic50, list("combination" = "before"))
args_ic50_after <- c(args_ic50, list("combination" = "after"))
args_ic80_before <- c(args_ic80, list("combination" = "before"))
args_ic80_after <- c(args_ic80, list("combination" = "after"))
args_sens_ic50_before <- c(args_sens_ic50, list("combination" = "before"))
args_sens_ic50_after <- c(args_sens_ic50, list("combination" = "after"))
args_sens_ic80_before <- c(args_sens_ic80, list("combination" = "before"))
args_sens_ic80_after <- c(args_sens_ic80, list("combination" = "after"))

set.seed(20230715)
cv_folds_continuous <- vimp::make_folds(y = analysis_dataset$ic50_additive, V = args$K)
# make stratified CV folds for combine-before-predict sens and individual-nAb sens
# need to do 5-fold CV for binary outcomes, since I have fewer than 10 in minority class
set.seed(20230719)
all_sens_names <- names(analysis_dataset %>% select(starts_with("sens_")))
for (i in 1:length(all_sens_names)) {
  this_y <- analysis_dataset %>% pull(all_sens_names[i])
  init_folds <- vimp::make_folds(y = this_y[!is.na(this_y)], V = binary_K, stratified = TRUE)
  na_folds <- vimp::make_folds(y = rep(1, sum(is.na(this_y))), V = binary_K)
  full_folds <- vector("numeric", length = nrow(analysis_dataset))
  full_folds[!is.na(this_y)] <- init_folds
  full_folds[is.na(this_y)] <- na_folds
  eval(parse(text = paste0(
    "cv_folds_", all_sens_names[i], " <- full_folds"
  )))
}

# set up a dataset that contains outcomes, predictions, folds for each observation (wide)
wide_pred_dataset <- analysis_dataset %>% 
  select(starts_with("ic"), starts_with("sens"))
for (i in 1:length(nice_bnabs)) {
  eval(parse(text = paste0("wide_pred_dataset$folds_sens_ic50_", nice_bnabs[i], " <- cv_folds_sens_ic50_", nice_bnabs[i])))
  eval(parse(text = paste0("wide_pred_dataset$folds_sens_ic80_", nice_bnabs[i], " <- cv_folds_sens_ic80_", nice_bnabs[i])))
  eval(parse(text = paste0("wide_pred_dataset$folds_ic50_", nice_bnabs[i], " <- cv_folds_continuous")))
  eval(parse(text = paste0("wide_pred_dataset$folds_ic80_", nice_bnabs[i], " <- cv_folds_continuous")))
}
wide_pred_dataset$folds_ic50_additive <- wide_pred_dataset$folds_ic50_bh <- cv_folds_continuous
wide_pred_dataset$folds_ic80_additive <- wide_pred_dataset$folds_ic80_bh <- cv_folds_continuous
wide_pred_dataset$folds_sens_ic50_additive <- cv_folds_sens_ic50_additive
wide_pred_dataset$folds_sens_ic50_bh <- cv_folds_sens_ic50_bh
wide_pred_dataset$folds_sens_ic80_additive <- cv_folds_sens_ic80_additive
wide_pred_dataset$folds_sens_ic80_bh <- cv_folds_sens_ic80_bh

long_folds_dataset <- wide_pred_dataset %>% 
  select(starts_with("folds")) %>% 
  mutate(id = row_number()) %>% 
  pivot_longer(cols = starts_with("folds"),
               names_to = c("outcome", "bnab"),
               names_pattern = "(folds_ic50|folds_ic80|folds_sens_ic50|folds_sens_ic80)_(.*)",
               values_to = "folds") %>% 
  mutate(outcome = gsub("folds_", "", outcome))

long_pred_dataset <- wide_pred_dataset %>% 
  select(-starts_with("folds")) %>% 
  mutate(id = row_number()) %>% 
  pivot_longer(cols = starts_with("ic") | starts_with("sens"), 
               names_to = c("outcome", "bnab"), names_pattern = "(ic50|ic80|sens_ic50|sens_ic80)_(.*)",
               values_to = "y") %>% 
  left_join(long_folds_dataset, by = c("id", "outcome", "bnab")) %>% 
  mutate(combo_time = case_when(
    bnab == "additive" | bnab == "bh" ~ "before",
    !(bnab == "additive" | bnab == "bh") ~ "after"
  )) %>% 
  mutate(method = case_when(bnab == "additive" ~ "additive",
                            bnab == "bh" ~ "Bliss-Hill",
                            !(bnab == "additive" | bnab == "bh") ~ "individual"),
         bnab = ifelse(bnab == "additive" | bnab == "bh", nice_bnab_combo, bnab)) %>% 
  select(-id) %>% 
  filter(!(method == "individual" & bnab == nice_bnab_combo))
long_pred_dataset$preds <- NA
# set up parallelization
num_cores <- parallel::detectCores()
cl <- parallel::makeCluster(num_cores)
doParallel::registerDoParallel(cl)

# get CV predictions of IC-50 --------------------------------------------------
# only have to run this once, since I've already combined
set.seed(20230713)
args_ic50_additive <- c(args_ic50_before, list("combination_method" = "additive"))
this_boolean <- long_pred_dataset$outcome == "ic50" & long_pred_dataset$bnab == nice_bnab_combo &
  long_pred_dataset$combo_time == "before" & long_pred_dataset$method == "additive"
these_folds <- long_pred_dataset$folds[this_boolean]
preds_ic50_additive <- get_one_bnab_preds(args_ic50_additive, analysis_dataset, bnabs[1],
                                          folds = these_folds,
                                          parallel = args$parallel)
long_pred_dataset[this_boolean, ]$preds[preds_ic50_additive$cc_indx] <- preds_ic50_additive$preds

set.seed(20230713)
args_ic50_bh <- c(args_ic50_before, list("combination_method" = "Bliss-Hill"))
this_boolean <- long_pred_dataset$outcome == "ic50" & long_pred_dataset$bnab == nice_bnab_combo &
  long_pred_dataset$combo_time == "before" & long_pred_dataset$method == "Bliss-Hill"
these_folds <- long_pred_dataset$folds[this_boolean]
preds_ic50_bh <- get_one_bnab_preds(args_ic50_bh, analysis_dataset, bnabs[1],
                                    folds = these_folds,
                                    parallel = args$parallel)
long_pred_dataset[this_boolean, ]$preds[preds_ic50_bh$cc_indx] <- preds_ic50_bh$preds

# need to run for each bnAb, then combine
preds_ic50 <- lapply(as.list(seq_len(length(bnabs))), function(j) {
  set.seed(20230713)
  this_boolean <- long_pred_dataset$outcome == "ic50" & long_pred_dataset$bnab == nice_bnabs[j] &
    long_pred_dataset$combo_time == "after"
  these_folds <- long_pred_dataset$folds[this_boolean]
  these_preds <- get_one_bnab_preds(args_ic50_after, analysis_dataset, bnabs[j],
                                    folds = these_folds,
                                    parallel = args$parallel)
  these_preds$bnab <- nice_bnabs[j]
  these_preds$method <- "individual"
  these_preds$combo_time <- "after"
  these_preds
})
preds_ic50_filled <- fill_preds_mat(preds = preds_ic50, n = nrow(analysis_dataset), 
                                    cv_folds = cv_folds_continuous)
for (i in seq_len(length(bnabs))) {
  this_boolean <- long_pred_dataset$outcome == "ic50" & long_pred_dataset$bnab == nice_bnabs[i] &
    long_pred_dataset$combo_time == "after" & long_pred_dataset$method == "individual"
  long_pred_dataset[this_boolean, ]$preds[preds_ic50[[i]]$cc_indx] <- preds_ic50[[i]]$preds
}
preds_mat_ic50 <- 10 ^ do.call(cbind, lapply(preds_ic50_filled, function(x) x$preds))

# get CV predictions of IC-80 --------------------------------------------------
set.seed(20230715)
args_ic80_additive <- c(args_ic80_before, list("combination_method" = "additive"))
this_boolean <- long_pred_dataset$outcome == "ic80" & long_pred_dataset$bnab == nice_bnab_combo &
  long_pred_dataset$combo_time == "before" & long_pred_dataset$method == "additive"
these_folds <- long_pred_dataset$folds[this_boolean]
preds_ic80_additive <- get_one_bnab_preds(args_ic80_additive, analysis_dataset, bnabs[1],
                                          folds = these_folds,
                                          parallel = args$parallel)
long_pred_dataset[this_boolean, ]$preds[preds_ic80_additive$cc_indx] <- preds_ic80_additive$preds

set.seed(20230715)
args_ic80_bh <- c(args_ic80_before, list("combination_method" = "Bliss-Hill"))
this_boolean <- long_pred_dataset$outcome == "ic80" & long_pred_dataset$bnab == nice_bnab_combo &
  long_pred_dataset$combo_time == "before" & long_pred_dataset$method == "Bliss-Hill"
these_folds <- long_pred_dataset$folds[this_boolean]
preds_ic80_bh <- get_one_bnab_preds(args_ic80_bh, analysis_dataset, bnabs[1],
                                    folds = these_folds,
                                    parallel = args$parallel)
long_pred_dataset[this_boolean, ]$preds[preds_ic80_bh$cc_indx] <- preds_ic80_bh$preds

preds_ic80 <- lapply(as.list(seq_len(length(bnabs))), function(j) {
  set.seed(20230715)
  this_boolean <- long_pred_dataset$outcome == "ic80" & long_pred_dataset$bnab == nice_bnabs[j] &
    long_pred_dataset$combo_time == "after"
  these_folds <- long_pred_dataset$folds[this_boolean]
  these_preds <- get_one_bnab_preds(args_ic80_after, analysis_dataset, bnabs[j],
                                    folds = these_folds,
                                    parallel = args$parallel)
  these_preds$bnab <- nice_bnabs[j]
  these_preds$method <- "individual"
  these_preds$combo_time <- "after"
  these_preds
})
preds_ic80_filled <- fill_preds_mat(preds = preds_ic80, n = nrow(analysis_dataset),
                                    cv_folds = cv_folds_continuous)
for (i in seq_len(length(bnabs))) {
  this_boolean <- long_pred_dataset$outcome == "ic80" & long_pred_dataset$bnab == nice_bnabs[i] &
    long_pred_dataset$combo_time == "after" & long_pred_dataset$method == "individual"
  long_pred_dataset[this_boolean, ]$preds[preds_ic80[[i]]$cc_indx] <- preds_ic80[[i]]$preds
}
preds_mat_ic80 <- 10 ^ do.call(cbind, lapply(preds_ic80_filled, function(x) x$preds))

# get CV predictions of susceptibility (IC50 < 1) ---------------------------------
set.seed(20230719)
args_sens_ic50_additive <- c(args_sens_ic50_before, list("combination_method" = "additive"))
this_boolean <- long_pred_dataset$outcome == "sens_ic50" & long_pred_dataset$bnab == nice_bnab_combo &
  long_pred_dataset$combo_time == "before" & long_pred_dataset$method == "additive"
these_folds <- long_pred_dataset$folds[this_boolean]
preds_sens_ic50_additive <- get_one_bnab_preds(args_sens_ic50_additive,
                                               analysis_dataset, bnabs[1],
                                               folds = these_folds,
                                               parallel = args$parallel)
long_pred_dataset[this_boolean, ]$preds[preds_sens_ic50_additive$cc_indx] <- preds_sens_ic50_additive$preds

set.seed(20230719)
args_sens_ic50_bh <- c(args_sens_ic50_before, list("combination_method" = "Bliss-Hill"))
this_boolean <- long_pred_dataset$outcome == "sens_ic50" & long_pred_dataset$bnab == nice_bnab_combo &
  long_pred_dataset$combo_time == "before" & long_pred_dataset$method == "Bliss-Hill"
these_folds <- long_pred_dataset$folds[this_boolean]
preds_sens_ic50_bh <- get_one_bnab_preds(args_sens_ic50_bh,
                                               analysis_dataset, bnabs[1],
                                               folds = these_folds,
                                               parallel = args$parallel)
long_pred_dataset[this_boolean, ]$preds[preds_sens_ic50_bh$cc_indx] <- preds_sens_ic50_bh$preds

# get individual bnAb predictions of susceptibility
preds_sens_ic50 <- lapply(as.list(seq_len(length(bnabs))), function(j) {
  set.seed(20230719)
  these_args <- args_sens_ic50_after
  these_args$outcome <- paste0(args_sens_ic50_after$outcome, "_", nice_bnabs[j])
  this_boolean <- long_pred_dataset$outcome == "sens_ic50" & long_pred_dataset$bnab == nice_bnabs[j] &
    long_pred_dataset$combo_time == "after"
  these_folds <- long_pred_dataset$folds[this_boolean]
  these_preds <- get_one_bnab_preds(args_sens_ic50_after, analysis_dataset, bnabs[j],
                                    folds = these_folds,
                                    parallel = args$parallel)
  these_preds$bnab <- nice_bnabs[j]
  these_preds$method <- "individual"
  these_preds$combo_time <- "after"
  these_preds
})
# preds_sens_ic50_filled <- fill_preds_mat(preds = preds_sens_ic50, n = nrow(analysis_dataset),
#                                          cv_folds = cv_folds_binary)
for (i in seq_len(length(bnabs))) {
  this_boolean <- long_pred_dataset$outcome == "sens_ic50" & long_pred_dataset$bnab == nice_bnabs[i] &
    long_pred_dataset$combo_time == "after" & long_pred_dataset$method == "individual"
  long_pred_dataset[this_boolean, ]$preds[preds_sens_ic50[[i]]$cc_indx] <- preds_sens_ic50[[i]]$preds
}

# get CV predictions of susceptibility (IC80 < 1) ---------------------------------
set.seed(20230719)
args_sens_ic80_additive <- c(args_sens_ic80_before, list("combination_method" = "additive"))
this_boolean <- long_pred_dataset$outcome == "sens_ic80" & long_pred_dataset$bnab == nice_bnab_combo &
  long_pred_dataset$combo_time == "before" & long_pred_dataset$method == "additive"
these_folds <- long_pred_dataset$folds[this_boolean]
preds_sens_ic80_additive <- get_one_bnab_preds(args_sens_ic80_additive,
                                               analysis_dataset, bnabs[1],
                                               folds = these_folds,
                                               parallel = args$parallel)
long_pred_dataset[this_boolean, ]$preds[preds_sens_ic80_additive$cc_indx] <- preds_sens_ic80_additive$preds

set.seed(20230719)
args_sens_ic80_bh <- c(args_sens_ic80_before, list("combination_method" = "Bliss-Hill"))
this_boolean <- long_pred_dataset$outcome == "sens_ic80" & long_pred_dataset$bnab == nice_bnab_combo &
  long_pred_dataset$combo_time == "before" & long_pred_dataset$method == "Bliss-Hill"
these_folds <- long_pred_dataset$folds[this_boolean]
preds_sens_ic80_bh <- get_one_bnab_preds(args_sens_ic80_bh,
                                         analysis_dataset, bnabs[1],
                                         folds = these_folds,
                                         parallel = args$parallel)
long_pred_dataset[this_boolean, ]$preds[preds_sens_ic80_bh$cc_indx] <- preds_sens_ic80_bh$preds

# get individual bnAb predictions of susceptibility
preds_sens_ic80 <- lapply(as.list(seq_len(length(bnabs))), function(j) {
  set.seed(20230719)
  these_args <- args_sens_ic80_after
  these_args$outcome <- paste0(args_sens_ic80_after$outcome, "_", nice_bnabs[j])
  this_boolean <- long_pred_dataset$outcome == "sens_ic80" & long_pred_dataset$bnab == nice_bnabs[j] &
    long_pred_dataset$combo_time == "after"
  these_folds <- long_pred_dataset$folds[this_boolean]
  these_preds <- get_one_bnab_preds(args_sens_ic80_after, analysis_dataset, bnabs[j],
                                    folds = these_folds,
                                    parallel = args$parallel)
  these_preds$bnab <- nice_bnabs[j]
  these_preds$method <- "individual"
  these_preds$combo_time <- "after"
  these_preds
})
# preds_sens_ic80_filled <- fill_preds_mat(preds = preds_sens_ic80, n = nrow(analysis_dataset),
#                                          cv_folds = cv_folds_binary)
for (i in seq_len(length(bnabs))) {
  this_boolean <- long_pred_dataset$outcome == "sens_ic80" & long_pred_dataset$bnab == nice_bnabs[i] &
    long_pred_dataset$combo_time == "after" & long_pred_dataset$method == "individual"
  long_pred_dataset[this_boolean, ]$preds[preds_sens_ic80[[i]]$cc_indx] <- preds_sens_ic80[[i]]$preds
}

parallel::stopCluster(cl)
# get combined CV predictions --------------------------------------------------
post_pred_ic50_additive <- log10(apply(preds_mat_ic50, 1, function(x) additive_combo_method(x)))
post_pred_ic80_additive <- log10(apply(preds_mat_ic80, 1, function(x) additive_combo_method(x)))
post_pred_ic50_bh <- log10(predict_all_bh(conc = 0.5, ic50 = preds_mat_ic50, 
                                          ic80 = preds_mat_ic80))
post_pred_ic80_bh <- log10(predict_all_bh(conc = 0.8, ic50 = preds_mat_ic50, 
                                          ic80 = preds_mat_ic80))
# predicted susceptibility based on predictions of individual IC50, IC80
post_sens_ic50_additive <- as.numeric(post_pred_ic50_additive < 0)
post_sens_ic80_additive <- as.numeric(post_pred_ic50_additive < 0)
post_sens_ic50_bh <- as.numeric(post_pred_ic50_bh < 0)
post_sens_ic80_bh <- as.numeric(post_pred_ic80_bh < 0)

# add to the dataset
not_na_ic50_additive <- !is.na(post_pred_ic50_additive)
post_pred_ic50_additive_df <- data.frame(
  "preds" = post_pred_ic50_additive[not_na_ic50_additive], 
  "folds" = cv_folds_continuous[not_na_ic50_additive],
  "cc_indx" = which(not_na_ic50_additive), "y" = analysis_dataset$ic50_additive[not_na_ic50_additive], 
  "outcome_type" = "ic50", "bnab" = nice_bnab_combo, "method" = "additive", "combo_time" = "after"
)
post_pred_sens_ic50_additive_df <- data.frame(
  "preds" = post_sens_ic50_additive[not_na_ic50_additive], 
  "folds" = cv_folds_continuous[not_na_ic50_additive],
  "cc_indx" = which(not_na_ic50_additive), "y" = analysis_dataset$sens_ic50_additive[not_na_ic50_additive], 
  "outcome_type" = "sens_ic50", "bnab" = nice_bnab_combo, "method" = "additive", "combo_time" = "after"
)

not_na_ic50_bh <- !is.na(post_pred_ic50_bh)
post_pred_ic50_bh_df <- data.frame(
  "preds" = post_pred_ic50_bh[not_na_ic50_bh], 
  "folds" = cv_folds_continuous[not_na_ic50_bh],
  "cc_indx" = which(not_na_ic50_bh), "y" = analysis_dataset$ic50_bh[not_na_ic50_bh], 
  "outcome_type" = "ic50", "bnab" = nice_bnab_combo, "method" = "Bliss-Hill", "combo_time" = "after"
)
post_pred_sens_ic50_bh_df <- data.frame(
  "preds" = post_sens_ic50_bh[not_na_ic50_bh], 
  "folds" = cv_folds_continuous[not_na_ic50_bh],
  "cc_indx" = which(not_na_ic50_bh), "y" = analysis_dataset$sens_ic50_bh[not_na_ic50_bh], 
  "outcome_type" = "sens_ic50", "bnab" = nice_bnab_combo, "method" = "Bliss-Hill", "combo_time" = "after"
)

post_pred_ic80_additive_df <- data.frame(
  "preds" = post_pred_ic80_additive[preds_ic80_additive$cc_indx], "folds" = preds_ic80_additive$folds,
  "cc_indx" = preds_ic80_additive$cc_indx, "y" = preds_ic80_additive$y, "outcome_type" = "ic80",
  "bnab" = nice_bnab_combo, "method" = "additive", "combo_time" = "after"
)
post_pred_sens_ic80_additive_df <- data.frame(
  "preds" = post_sens_ic80_additive[preds_ic80_additive$cc_indx], 
  "folds" = cv_folds_continuous[preds_ic80_additive$cc_indx],
  "cc_indx" = preds_ic80_additive$cc_indx, "y" = analysis_dataset$sens_ic80_additive[preds_ic80_additive$cc_indx], 
  "outcome_type" = "sens_ic80", "bnab" = nice_bnab_combo, "method" = "additive", "combo_time" = "after"
)
post_pred_ic80_bh_df <- data.frame(
  "preds" = post_pred_ic80_bh[preds_ic80_bh$cc_indx], "folds" = preds_ic80_bh$folds,
  "cc_indx" = preds_ic80_bh$cc_indx, "y" = preds_ic80_bh$y, "outcome_type" = "ic80",
  "bnab" = nice_bnab_combo, "method" = "Bliss-Hill", "combo_time" = "after"
)
post_pred_sens_ic80_bh_df <- data.frame(
  "preds" = post_sens_ic80_bh[preds_ic80_bh$cc_indx], 
  "folds" = cv_folds_continuous[preds_ic80_bh$cc_indx],
  "cc_indx" = preds_ic80_bh$cc_indx, "y" = analysis_dataset$sens_ic80_bh[preds_ic80_bh$cc_indx], 
  "outcome_type" = "sens_ic80", "bnab" = nice_bnab_combo, "method" = "Bliss-Hill", "combo_time" = "after"
)

all_post_preds <- do.call(rbind, args = list(
  post_pred_ic50_additive_df, post_pred_ic50_bh_df,
  post_pred_ic80_additive_df, post_pred_ic80_bh_df,
  post_pred_sens_ic50_additive_df, post_pred_sens_ic50_bh_df,
  post_pred_sens_ic80_additive_df, post_pred_sens_ic80_bh_df
)) %>% 
  rename(outcome = outcome_type) %>% 
  select(outcome, bnab, y, combo_time, method, folds, preds)


# if we have lab-measured combo neutralization, want prediction performance
# against this too
if (has_combo_neut) {
  combo_ic50 <- analysis_dataset %>% 
    pull(grep(paste0("ic50_", nice_bnab_combo), names(analysis_dataset))[1])
  combo_ic80 <- analysis_dataset %>% 
    pull(grep(paste0("ic80_", nice_bnab_combo), names(analysis_dataset))[1])
  combo_sens_ic50 <- analysis_dataset %>% 
    pull(grep(paste0("sens_ic50_", nice_bnab_combo), names(analysis_dataset)))
  combo_sens_ic80 <- analysis_dataset %>% 
    pull(grep(paste0("sens_ic80_", nice_bnab_combo), names(analysis_dataset)))
  # additive model, pre-prediction combination
  this_boolean <- long_pred_dataset$outcome == "ic50" & long_pred_dataset$bnab == nice_bnab_combo &
    long_pred_dataset$method == "additive"
  these_preds <- long_pred_dataset[this_boolean, ]
  preds_ic50_additive_true_y <- data.frame(
    "preds" = these_preds$preds[!is.na(combo_ic50)],
    "folds" = these_preds$folds[!is.na(combo_ic50)],
    "cc_indx" = which(!is.na(combo_ic50)), "y" = combo_ic50[!is.na(combo_ic50)],
    "outcome_type" = "true_combo_ic50", bnab = nice_bnab_combo, "method" = "additive",
    "combo_time" = "before"
  )
  this_boolean <- long_pred_dataset$outcome == "ic80" & long_pred_dataset$bnab == nice_bnab_combo &
    long_pred_dataset$method == "additive"
  these_preds <- long_pred_dataset[this_boolean, ]
  preds_ic80_additive_true_y <- data.frame(
    "preds" = these_preds$preds[!is.na(combo_ic80)],
    "folds" = these_preds$folds[!is.na(combo_ic80)],
    "cc_indx" = which(!is.na(combo_ic80)), "y" = combo_ic80[!is.na(combo_ic80)],
    "outcome_type" = "true_combo_ic80", bnab = nice_bnab_combo, "method" = "additive",
    "combo_time" = "before"
  )
  this_boolean <- long_pred_dataset$outcome == "sens_ic50" & long_pred_dataset$bnab == nice_bnab_combo &
    long_pred_dataset$method == "additive"
  these_preds <- long_pred_dataset[this_boolean, ]
  preds_sens_ic50_additive_true_y <- data.frame(
    "preds" = these_preds$preds[!is.na(combo_sens_ic50)],
    "folds" = these_preds$folds[!is.na(combo_sens_ic50)],
    "cc_indx" = which(!is.na(combo_sens_ic50)), "y" = combo_sens_ic50[!is.na(combo_sens_ic50)],
    "outcome_type" = "true_combo_sens_ic50", bnab = nice_bnab_combo, "method" = "additive",
    "combo_time" = "before"
  )
  this_boolean <- long_pred_dataset$outcome == "sens_ic80" & long_pred_dataset$bnab == nice_bnab_combo &
    long_pred_dataset$method == "additive"
  these_preds <- long_pred_dataset[this_boolean, ]
  preds_sens_ic80_additive_true_y <- data.frame(
    "preds" = these_preds$preds[!is.na(combo_sens_ic80)],
    "folds" = these_preds$folds[!is.na(combo_sens_ic80)],
    "cc_indx" = which(!is.na(combo_sens_ic80)), "y" = combo_sens_ic80[!is.na(combo_sens_ic80)],
    "outcome_type" = "true_combo_sens_ic80", bnab = nice_bnab_combo, "method" = "additive",
    "combo_time" = "before"
  )  
  # BH model, pre-prediction combination
  this_boolean <- long_pred_dataset$outcome == "ic50" & long_pred_dataset$bnab == nice_bnab_combo &
    long_pred_dataset$method == "Bliss-Hill"
  these_preds <- long_pred_dataset[this_boolean, ]
  preds_ic50_bh_true_y <- data.frame(
    "preds" = these_preds$preds[!is.na(combo_ic50)],
    "folds" = these_preds$folds[!is.na(combo_ic50)],
    "cc_indx" = which(!is.na(combo_ic50)), "y" = combo_ic50[!is.na(combo_ic50)],
    "outcome_type" = "true_combo_ic50", bnab = nice_bnab_combo, "method" = "Bliss-Hill",
    "combo_time" = "before"
  )
  this_boolean <- long_pred_dataset$outcome == "ic80" & long_pred_dataset$bnab == nice_bnab_combo &
    long_pred_dataset$method == "Bliss-Hill"
  these_preds <- long_pred_dataset[this_boolean, ]
  preds_ic80_bh_true_y <- data.frame(
    "preds" = these_preds$preds[!is.na(combo_ic80)],
    "folds" = these_preds$folds[!is.na(combo_ic80)],
    "cc_indx" = which(!is.na(combo_ic80)), "y" = combo_ic80[!is.na(combo_ic80)],
    "outcome_type" = "true_combo_ic80", bnab = nice_bnab_combo, "method" = "Bliss-Hill",
    "combo_time" = "before"
  )
  this_boolean <- long_pred_dataset$outcome == "sens_ic50" & long_pred_dataset$bnab == nice_bnab_combo &
    long_pred_dataset$method == "Bliss-Hill"
  these_preds <- long_pred_dataset[this_boolean, ]
  preds_sens_ic50_bh_true_y <- data.frame(
    "preds" = these_preds$preds[!is.na(combo_sens_ic50)],
    "folds" = these_preds$folds[!is.na(combo_sens_ic50)],
    "cc_indx" = which(!is.na(combo_sens_ic50)), "y" = combo_sens_ic50[!is.na(combo_sens_ic50)],
    "outcome_type" = "true_combo_sens_ic50", bnab = nice_bnab_combo, "method" = "Bliss-Hill",
    "combo_time" = "before"
  )
  this_boolean <- long_pred_dataset$outcome == "sens_ic80" & long_pred_dataset$bnab == nice_bnab_combo &
    long_pred_dataset$method == "Bliss-Hill"
  these_preds <- long_pred_dataset[this_boolean, ]
  preds_sens_ic80_bh_true_y <- data.frame(
    "preds" = these_preds$preds[!is.na(combo_sens_ic80)],
    "folds" = these_preds$folds[!is.na(combo_sens_ic80)],
    "cc_indx" = which(!is.na(combo_sens_ic80)), "y" = combo_sens_ic80[!is.na(combo_sens_ic80)],
    "outcome_type" = "true_combo_sens_ic80", bnab = nice_bnab_combo, "method" = "Bliss-Hill",
    "combo_time" = "before"
  )  
  # additive model, post-prediction combination
  post_preds_ic50_additive_true_y <- data.frame(
    "preds" = post_pred_ic50_additive[!is.na(combo_ic50)],
    "folds" = cv_folds_continuous[!is.na(combo_ic50)],
    "cc_indx" = which(!is.na(combo_ic50)), "y" = combo_ic50[!is.na(combo_ic50)],
    "outcome_type" = "true_combo_ic50", bnab = nice_bnab_combo, "method" = "additive",
    "combo_time" = "after"
  )
  post_preds_ic80_additive_true_y <- data.frame(
    "preds" = post_pred_ic80_additive[!is.na(combo_ic80)],
    "folds" = cv_folds_continuous[!is.na(combo_ic80)],
    "cc_indx" = which(!is.na(combo_ic80)), "y" = combo_ic80[!is.na(combo_ic80)],
    "outcome_type" = "true_combo_ic80", bnab = nice_bnab_combo, "method" = "additive",
    "combo_time" = "after"
  )
  post_preds_sens_ic50_additive_true_y <- data.frame(
    "preds" = post_sens_ic50_additive[!is.na(combo_sens_ic50)],
    "folds" = cv_folds_continuous[!is.na(combo_sens_ic50)],
    "cc_indx" = which(!is.na(combo_sens_ic50)), "y" = combo_sens_ic50[!is.na(combo_sens_ic50)],
    "outcome_type" = "true_combo_sens_ic50", bnab = nice_bnab_combo, "method" = "additive",
    "combo_time" = "after"
  )
  post_preds_sens_ic80_additive_true_y <- data.frame(
    "preds" = post_sens_ic80_additive[!is.na(combo_sens_ic80)],
    "folds" = cv_folds_continuous[!is.na(combo_sens_ic80)],
    "cc_indx" = which(!is.na(combo_sens_ic80)), "y" = combo_sens_ic80[!is.na(combo_sens_ic80)],
    "outcome_type" = "true_combo_sens_ic80", bnab = nice_bnab_combo, "method" = "additive",
    "combo_time" = "after"
  )  
  # BH model, post-prediction combination
  post_preds_ic50_bh_true_y <- data.frame(
    "preds" = post_pred_ic50_bh[!is.na(combo_ic50)],
    "folds" = cv_folds_continuous[!is.na(combo_ic50)],
    "cc_indx" = which(!is.na(combo_ic50)), "y" = combo_ic50[!is.na(combo_ic50)],
    "outcome_type" = "true_combo_ic50", bnab = nice_bnab_combo, "method" = "Bliss-Hill",
    "combo_time" = "after"
  )
  post_preds_ic80_bh_true_y <- data.frame(
    "preds" = post_pred_ic80_bh[!is.na(combo_ic80)],
    "folds" = cv_folds_continuous[!is.na(combo_ic80)],
    "cc_indx" = which(!is.na(combo_ic80)), "y" = combo_ic80[!is.na(combo_ic80)],
    "outcome_type" = "true_combo_ic80", bnab = nice_bnab_combo, "method" = "Bliss-Hill",
    "combo_time" = "after"
  )
  post_preds_sens_ic50_bh_true_y <- data.frame(
    "preds" = post_sens_ic50_bh[!is.na(combo_sens_ic50)],
    "folds" = cv_folds_continuous[!is.na(combo_sens_ic50)],
    "cc_indx" = which(!is.na(combo_sens_ic50)), "y" = combo_sens_ic50[!is.na(combo_sens_ic50)],
    "outcome_type" = "true_combo_sens_ic50", bnab = nice_bnab_combo, "method" = "Bliss-Hill",
    "combo_time" = "after"
  )
  post_preds_sens_ic80_bh_true_y <- data.frame(
    "preds" = post_sens_ic80_bh[!is.na(combo_sens_ic80)],
    "folds" = cv_folds_continuous[!is.na(combo_sens_ic80)],
    "cc_indx" = which(!is.na(combo_sens_ic80)), "y" = combo_sens_ic80[!is.na(combo_sens_ic80)],
    "outcome_type" = "true_combo_sens_ic80", bnab = nice_bnab_combo, "method" = "Bliss-Hill",
    "combo_time" = "after"
  )  
  all_true_preds <- rbind.data.frame(
    preds_ic50_additive_true_y, preds_ic80_additive_true_y, 
    preds_sens_ic50_additive_true_y, preds_sens_ic80_additive_true_y, 
    preds_ic50_bh_true_y, preds_ic80_bh_true_y, 
    preds_sens_ic50_bh_true_y, preds_sens_ic80_bh_true_y, 
    post_preds_ic50_additive_true_y, post_preds_ic80_additive_true_y, 
    post_preds_sens_ic50_additive_true_y, post_preds_sens_ic80_additive_true_y, 
    post_preds_ic50_bh_true_y, post_preds_ic80_bh_true_y, 
    post_preds_sens_ic50_bh_true_y, post_preds_sens_ic80_bh_true_y
  ) %>% 
    rename(outcome = outcome_type) %>% 
    select(outcome, bnab, y, combo_time, method, folds, preds)
} else {
  all_true_preds <- NULL
}

# get all predictions together, along with outcomes and CV folds; this is a long dataset
all_preds <- rbind(long_pred_dataset, all_post_preds, all_true_preds)

# save
if (!dir.exists(here::here("..", "results"))) {
  dir.create(here::here("..", "results"), recursive = TRUE)
}
saveRDS(all_preds, file = here::here("..", "results", paste0("preds_", nice_bnab_combo, ".rds")))
print("Predictions complete.")
