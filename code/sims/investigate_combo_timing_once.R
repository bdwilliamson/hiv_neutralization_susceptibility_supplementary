# one run of the investigation

investigate_combo_timing_once <- function(mc_id = 1, n = 100, n_mabs = 3, p = 1000, model = "simple", correlated = FALSE, strength = "strong") {
  # generate data for each antibody
  dat <- gen_data(n = n, p = p, correlated = correlated, n_mabs = n_mabs, strength = strength)
  x_mat <- as.matrix(dat$x)
  # get combined ic80, susceptibility using additive model
  all_ic80 <- 10 ^ do.call(cbind, dat$y_cont)
  ic80_additive <- log10(apply(all_ic80, 1, function(x) additive_combo_method(x)))
  susc_additive <- as.numeric(ic80_additive < 0)
  # fit a cross-validated lasso to each antibody's data individually
  continuous_folds <- vimp::make_folds(y = ic80_additive, V = 10, stratified = FALSE)
  binary_folds <- vimp::make_folds(y = susc_additive, V = 10, stratified = TRUE)
  ind_cont_cv_preds <- lapply(as.list(1:n_mabs), function(i) {
    get_cv_preds(y = dat$y_cont[, i], x = x_mat, K = 10, parallel = FALSE,
                 folds = continuous_folds)
  })
  ind_bin_cv_preds <- lapply(as.list(1:n_mabs), function(i) {
    get_cv_preds(y = dat$y_bin[, i], x = x_mat, K = 10, parallel = FALSE,
                 folds = binary_folds)
  })
  # combine predictions
  all_cont_preds <- 10 ^ do.call(cbind, lapply(ind_cont_cv_preds, function(pred) pred$preds))
  post_pred_comb <- log10(apply(all_cont_preds, 1, function(x) additive_combo_method(x)))
  post_pred_comb_sens <- as.numeric(post_pred_comb < 0)
  # fit cv lasso to combination values
  comb_cont_cv_preds <- get_cv_preds(y = ic80_additive, x = x_mat,
                                     K = 10, parallel = FALSE,
                                     folds = continuous_folds)
  comb_bin_cv_preds <- get_cv_preds(y = susc_additive, x = x_mat,
                                    K = 10, parallel = FALSE,
                                    folds = binary_folds)
  
  # get prediction performance
  ind_cont_perf <- lapply(as.list(1:n_mabs), function(i) {
    cvR2(predictions = ind_cont_cv_preds[[i]]$preds, y = dat$y_cont[, i],
         folds = continuous_folds)
  })
  pre_pred_comb_cont_perf <- cvR2(predictions = comb_cont_cv_preds$preds,
                                  y = ic80_additive, folds = continuous_folds)
  post_pred_comb_cont_perf <- cvR2(predictions = post_pred_comb, y = ic80_additive,
                                   folds = continuous_folds)
  ind_bin_perf <- lapply(as.list(1:n_mabs), function(i) {
    tryCatch(cvAUC::ci.cvAUC(predictions = ind_bin_cv_preds[[i]]$preds, labels = dat$y_bin[, i],
                             folds = binary_folds), error = function(e) {
           # remove any folds with zero 1s or zero 0s
           fold_nums <- unlist(lapply(as.list(sort(unique(binary_folds))), function(fold) {
             sum(dat$y_bin[binary_folds == fold, i])
           }))
           drop_folds <- which(fold_nums == 0)
           bool <- !(binary_folds %in% drop_folds)
           cvAUC::ci.cvAUC(predictions = ind_bin_cv_preds[[i]]$preds[bool],
                           labels = dat$y_bin[, i][bool],
                           folds = binary_folds[bool])
         })
  })
  pre_pred_comb_bin_perf <- cvAUC::ci.cvAUC(predictions = comb_bin_cv_preds$preds, 
                                     labels = susc_additive, folds = binary_folds)
  post_pred_comb_bin_perf <- cvAUC::ci.cvAUC(predictions = post_pred_comb_sens, labels = susc_additive,
                                      folds = binary_folds)
  # create output tibble
  point_ests_cont <- c(do.call(c, lapply(ind_cont_perf, function(x) x$cvR2)), pre_pred_comb_cont_perf$cvR2,
                       post_pred_comb_cont_perf$cvR2)
  point_ests_bin <- c(do.call(c, lapply(ind_bin_perf, function(x) x$cvAUC)), pre_pred_comb_bin_perf$cvAUC,
                      post_pred_comb_bin_perf$cvAUC)
  cis_cont <- rbind(do.call(rbind, lapply(ind_cont_perf, function(x) x$ci)),
                    pre_pred_comb_cont_perf$ci, post_pred_comb_cont_perf$ci)
  cis_bin <- rbind(do.call(rbind, lapply(ind_bin_perf, function(x) x$ci)),
                   pre_pred_comb_bin_perf$ci, post_pred_comb_bin_perf$ci)
  output <- tibble::tibble(mc_id = mc_id, n = n, n_mabs = n_mabs, p = p, 
                           model = model, correlated = correlated, strength = strength,
                           outcome = c(rep("IC80", length(point_ests_cont)),
                                       rep("Susceptibility", length(point_ests_bin))),
                           measure = c(rep("cvR2", length(point_ests_cont)),
                                       rep("cvAUC", length(point_ests_bin))),
                           mab = rep(c(paste0("ind_", 1:(length(point_ests_cont) - 2)), 
                                    "comb_pre", "comb_post"), 2),
                           point_est = c(point_ests_cont, point_ests_bin),
                           cil = c(cis_cont[, 1], cis_bin[, 1]),
                           ciu = c(cis_cont[, 2], cis_bin[, 2]))
  # return
  return(output)
}