# useful functions

# for combining neutralization data --------------------------------------------

get_iip <- function(ic50, ic80) {
  m <- log10(4) / (log10(ic80) - log10(ic50))
  return(10 ^ m / (ic50 ^ m + 10 ^ m))
}
# Additive method of Wagh et al. (2016)
# @param x a matrix of neutralization values
# @return a vector of combined neutralization values following the additive method of Wagh et al. (2016)
additive_combo_method <- function(x) {
  return(1 / sum(1 / x))
}

# compute individual mAb neutralization curve based on concentration and IC-50
# @param conc the concentration
# @param ic50 the IC-50 value
# @param m the slope (in our analyses, taken to be log(4) / [log(IC-80) - log(IC-50)])
individual_mab_neutralization <- function(conc, mab_conc, ic50, m) {
  f_mab_c <- (conc ^ m) / (ic50 ^ m + mab_conc ^ m)
  f_mab_c
}
# compute neutralization curves using the Bliss-Hill model (with independence assumption)
# @param conc the concentration
# @param ic50 vector or matrix of IC-50 values (columns are mAbs)
# @param ic80 vector or matrix of IC-80 values (columns are mAbs)
bliss_hill_predictions <- function(conc, ic50, ic80) {
  m <- log10(4) / (log10(ic80) - log10(ic50))
  if (!is.null(ncol(ic50))) {
    num_indices <- ncol(ic50)
    this_function <- function(indx, conc) {
      individual_mab_neutralization(conc = conc, mab_conc = conc / ncol(ic50), ic50[, indx], m[, indx])
    }
  } else {
    num_indices <- length(ic50)
    this_function <- function(indx, conc) {
      individual_mab_neutralization(conc = conc, mab_conc = conc / length(ic50), ic50[indx], m[indx])
    }
  }
  one_minus_fs <- do.call(cbind, sapply(1:num_indices, function(indx) 1 - this_function(indx, conc),
                                        simplify = FALSE))
  f_c <- 1 - apply(one_minus_fs, 1, prod)
  f_c
}

# Bliss-Hill method of Wagh et al. (2016)
# back-solve to obtain combination IC-50 or IC-80 from the Bliss-Hill model (with independence assumption)
# @param conc the concentration of interest
# @param ic50 the IC-50 values for all mAbs of interest (a vector)
# @param ic80 the IC-80 values for all mAbs of interest (a vector)
# @return the predicted IC (IC-50 if conc = 0.5, IC-80 if conc = 0.8)
predict_bh_concentration <- function(conc = 0.5, ic50, ic80) {
  optim_func <- function(another_conc, ic50, ic80, conc) {
    abs(bliss_hill_predictions(another_conc, ic50, ic80) - conc)
  }
  if (conc == 0.5) {
    init_par <- 1 / sum(1 / ic50)
    upper_lim <- switch(all(is.na(ic50)) + 1, min(ic50, na.rm = TRUE), NA)
  } else if (conc == 0.8) {
    init_par <- 1 / sum(1 / ic80)
    upper_lim <- switch(all(is.na(ic80)) + 1, min(ic80, na.rm = TRUE), NA)
  } else {
    init_par <- 0
    upper_lim <- 200
  }
  # if any individual one is NA, the whole thing should be
  if (any(is.na(c(ic50, ic80))) | is.na(upper_lim)) {
    ret <- NA
  } else { # do the optimization
    suppressWarnings(optimized <- optim(init_par, fn = optim_func, ic50 = ic50, ic80 = ic80, conc = conc,
                                        method = "Brent", lower = 0, upper = upper_lim,
                                        control = list(abstol = 1e-3, reltol = 1e-5)))
    ret <- optimized$par
  }
  ret
}
# Get Bliss-Hill predictions on a full dataset
# @param conc the concentration of interest
# @param ic50 a matrix of ic50 values
# @param ic80 a matrix of ic80 values
# @return the Bliss-Hill predictions for the entire matrix
predict_all_bh <- function(conc = 0.5, ic50 = NULL, ic80 = NULL) {
  sapply(1:nrow(ic50), function(i) {
    predict_bh_concentration(conc = conc, ic50 = ic50[i, ], ic80 = ic80[i, ])
  })
}

# fill NAs in the hxb2 columns with zeros
fill_na_hxb2 <- function(dt) {
  has_hxb2 <- names(dt)[grep("hxb2", names(dt))]
  for (j in has_hxb2) {
    dt[is.na(get(j)), (j):=0]
  }
  dt
}

# for getting predictions ------------------------------------------------------
# get cross-validated predictions for a given outcome
# @param y the outcome
# @param x the features
# @param K the number of cross-fitting folds
# @param folds the folds (if supplied)
# @return the cross-validated predictions
get_cv_preds <- function(y = NULL, x = NULL, K = 10, folds = NULL, parallel = TRUE) {
  is_binary <- length(unique(y)) == 2
  if (is_binary) {
    fam <- "binomial"
    this_nfolds <- K
  } else {
    fam <- "gaussian"
    this_nfolds <- 10
  }
  if (all(is.null(folds))) {
    folds <- vimp::make_folds(y = y, V = K, stratified = is_binary)  
  }
  cv_preds <- vector("numeric", length = length(y))
  for (k in seq_len(K)) {
    # break the data up
    train_y <- y[folds != k]
    train_x <- x[folds != k, ]
    test_y <- y[folds == k]
    test_x <- x[folds == k, ]
    # fit a lasso to the training data
    if (is_binary) {
      these_nfolds <- min(this_nfolds, min(table(train_y)))  
    } else {
      these_nfolds <- this_nfolds
    }
    if (these_nfolds < 3) {
      cv_folds_train <- vimp::make_folds(y = train_y, V = 10, stratified = is_binary)
      lasso_fit <- glmnet::cv.glmnet(x = train_x, y = train_y, family = "gaussian",
                                     nfolds = 10, parallel = parallel,
                                     foldid = cv_folds_train)
      
    } else {
      cv_folds_train <- vimp::make_folds(y = train_y, V = these_nfolds, stratified = is_binary)
      lasso_fit <- glmnet::cv.glmnet(x = train_x, y = train_y, family = fam,
                                     nfolds = these_nfolds, parallel = parallel,
                                     foldid = cv_folds_train)
    }
    # evaluate predictions on the testing data
    cv_preds[folds == k] <- predict(lasso_fit, newx = test_x,
                                    type = "response", s = "lambda.min")
  }
  return(list("preds" = cv_preds, "folds" = folds))
}

# get one set of predictions for a given bnab
get_one_bnab_preds <- function(args, analysis_dataset, bnab, folds = NULL, parallel = FALSE) {
  nice_nab_name <- janitor::make_clean_names(bnab)
  is_sens <- grepl("sens", args$outcome)
  if (args$combination == "before") {
    if (args$outcome == "ic50") {
      y <- switch(as.numeric(args$combination_method == "Bliss-Hill") + 1, analysis_dataset$ic50_additive, analysis_dataset$ic50_bh)
    } else if (args$outcome == "ic80") {
      y <- switch(as.numeric(args$combination_method == "Bliss-Hill") + 1, analysis_dataset$ic80_additive, analysis_dataset$ic80_bh)
    } else { # only can be "sens", since not doing multsens, since it only relies on individual nAbs
      if (grepl("ic50", args$outcome)) {
        y <- switch(as.numeric(args$combination_method == "Bliss-Hill") + 1, analysis_dataset$sens_ic50_additive, analysis_dataset$sens_ic50_bh)
      } else {
        y <- switch(as.numeric(args$combination_method == "Bliss-Hill") + 1, analysis_dataset$sens_ic80_additive, analysis_dataset$sens_ic80_bh)
      }
    }  
  } else {
    this_outcome <- which(grepl(args$outcome, names(analysis_dataset)) & grepl(nice_nab_name, names(analysis_dataset)))
    if (length(this_outcome) > 1 & !is_sens) {
      this_outcome <- which(grepl(args$outcome, names(analysis_dataset)) & grepl(nice_nab_name, names(analysis_dataset)) & !grepl("sens", names(analysis_dataset)))
    }
    if (length(this_outcome) > 1) {
      all_nms <- names(analysis_dataset)[this_outcome]
      bnabs_only <- gsub(paste0(args$outcome, "_"), "", all_nms)
      this_outcome <- this_outcome[startsWith(bnabs_only, nice_nab_name) & endsWith(bnabs_only, nice_nab_name)]
    }
    y <- analysis_dataset[[this_outcome]]
  }
  
  x <- as.matrix(
    analysis_dataset %>% 
      select(-starts_with("seq.id"), -starts_with("ic"), -starts_with("iip"), -starts_with("sens"))
  )
  cc <- complete.cases(y, x)
  y_cc <- y[cc]
  x_cc <- x[cc, ]
  folds_cc <- switch(as.numeric(length(folds) == length(y_cc)) + 1, folds[cc], folds)
  preds <- get_cv_preds(y = y_cc, x = x_cc, K = args$K, folds = folds_cc, parallel = parallel)
  return(list("preds" = preds$preds, "folds" = preds$folds, "cc_indx" = which(cc),
              "y" = y_cc, "outcome_type" = args$outcome))
}

# for filling a matrix of predictions to have the correct number of rows
# @param preds a list, where each element is the predictions for a given bnAb
# @param n the sample size for the entire analysis dataset
# @param cv_folds the full-dataset cross-validation folds
fill_preds_mat <- function(preds = NULL, n = 100, cv_folds = rep(1, n)) {
  preds_list <- vector("list", length = length(preds))
  for (j in 1:length(preds)) {
    preds_mat <- do.call(cbind.data.frame, preds[[j]])
    na_mat <- matrix(NA, nrow = n - nrow(preds_mat), ncol = ncol(preds_mat))
    colnames(na_mat) <- names(preds_mat)
    na_df <- as.data.frame(na_mat)
    if (nrow(na_df) != 0) {
      na_df$cc_indx <- (1:n)[-preds_mat$cc_indx]
      na_df$folds <- cv_folds[-preds_mat$cc_indx]
      na_df$bnab <- preds_mat$bnab[1]
      na_df$method <- preds_mat$method[1]
      na_df$outcome_type <- preds_mat$outcome_type[1]
      na_df$combo_time <- preds_mat$combo_time[1]
    }
    preds_mat_2 <- rbind.data.frame(preds_mat, na_df)
    preds_mat_2 <- preds_mat_2[order(preds_mat_2$cc_indx), ]
    preds_list[[j]] <- preds_mat_2
  }
  return(preds_list)
}

# for getting prediction performance -------------------------------------------
# CV-R^2
# @param predictions the predictions
# @param y the outcomes
# @param folds the folds
# @param ... other arguments to pass to vimp::measure_r_squared
# @return the cross-validated R-squared
cvR2 <- function(predictions = NULL, y = NULL, folds = NULL, ...) {
  outcome_var_lst <- vimp::measure_mse(fitted_values = rep(mean(y, na.rm = TRUE), length(y)), 
                                       y = y, ...)
  outcome_var <- outcome_var_lst$point_est
  outcome_var_eif <- outcome_var_lst$eif
  fold_mse_list <- lapply(as.list(1:length(unique(folds))), function(k) {
    vimp::measure_mse(fitted_values = predictions[folds == k],
                            y = y[folds == k], full_y = y, ...)
  })
  fold_mse <- unlist(lapply(fold_mse_list, function(fold_list) fold_list$point_est))
  fold_mse_eifs <- lapply(fold_mse_list, function(fold_list) fold_list$eif)
  fold_r2_eifs <- lapply(as.list(1:length(fold_mse_list)), function(k) {
    (-1) * as.vector(
      matrix(c(1 / outcome_var, (-1) * fold_mse[k] / (outcome_var ^ 2)), nrow = 1) %*%
        t(cbind(fold_mse_eifs[[k]], outcome_var_eif[folds == k]))
    )
  })
  fold_r2_var <- unlist(lapply(fold_r2_eifs, function(eif) mean(eif ^ 2, ...)))
  # fold_mse_var <- unlist(lapply(fold_mse_list, function(fold_list) mean(fold_list$eif ^ 2, ...)))
  cvr2 <- 1 - mean(fold_mse) / outcome_var
  cvr2_var <- mean(fold_r2_var)
  cvr2_ci <- vimp::vimp_ci(est = cvr2, se = sqrt(cvr2_var / length(y)), truncate = FALSE)
  return(list("fold_R2" = 1 - fold_mse / outcome_var, "cvR2" = cvr2, "ci" = cvr2_ci))
}

# functions for plotting and results -------------------------------------------
plot_nab_name <- function(bnab, bnabs, nice_bnabs) {
  final_bnab <- bnab
  for (i in 1:length(bnabs)) {
    final_bnab <- gsub(nice_bnabs[i], bnabs[i], final_bnab, fixed = TRUE)
  }
  return(final_bnab)
}