# generate data for the combination timing experiment

# @param n the sample size
# @param p the number of covariates
# @param correlated should the features be correlated?
# @return a dataset with n rows and p covariates (5 continuous, p - 5 binary)
gen_x <- function(n = 100, p = 1000, correlated = FALSE) {
  # generate several MVN covariates (like length, etc.) and many binary covariates
  if (!correlated) {
    Sigma <- diag(p)
  } else {
    Sigma2 <- 0.1 * (1 - diag(p))
    Sigma <- diag(p) + Sigma2
  }
  latent_x <- mvtnorm::rmvnorm(n, rep(0, p), sigma = Sigma)
  x <- latent_x 
  p_x <- 0.4
  x[, 6:p] <- apply(x[, 6:p], 2, function(z) as.numeric(z < quantile(z, p_x)))
  return(x)
}

# @param n the sample size
# @param x the covariates
# @param model "simple" or "complex" (for now, only simple)
# @param outcome_type "binary" or "continuous"
# @param n_mabs the number of mAbs to generate outcomes for
# @param strength the strength of association
gen_y <- function(n = 100, x = NULL, model = "simple", n_mabs = 3, strength = "strong") {
  p <- ncol(x)
  if (model != "simple") {
    # do nothing
  } else {
    if (strength != "weak") {
      beta1 <- rep(1, 10)
      beta2 <- rep(2, 10)
      beta3 <- rep(0.5, 10)
    } else {
      beta1 <- rep(0.2, 10)
      beta2 <- rep(0.4, 10)
      beta3 <- rep(0.1, 10)
    }
    beta_01 <- c(beta1, rep(0, p - 10))
    beta_02 <- c(beta2, rep(0, p - 10))
    beta_03 <- c(beta3, rep(0, p - 10))
    beta_0 <- list(beta_01, beta_02, beta_03)
  }
  ics <- vector("list", length = n_mabs)
  for (i in 1:n_mabs) {
    linear_predictor <- x %*% as.matrix(beta_0[[i]])
    # in this model, log10 ic80 ~ N(x\beta, 1)
    ics[[i]] <- linear_predictor + rnorm(n, 0, 1)
  }
  # since on log10 scale
  susc <- lapply(ics, function(z) as.numeric(z < 0))
  ic_df <- do.call(cbind.data.frame, ics)
  names(ic_df) <- paste0("ic_", 1:n_mabs)
  susc_df <- do.call(cbind.data.frame, susc)
  names(susc_df) <- paste0("susc_", 1:n_mabs)
  return(list("continuous" = ic_df, "binary" = susc_df))
}

# @param n the sample size
# @param p the number of covariates
# @param correlated should the features be correlated?
# @param model "simple" or "complex" (for now, only simple)
# @param outcome_type "binary" or "continuous"
# @param n_mabs the number of mAbs
# @param strength the outcome-feature strength
# @return a full dataset with n rows, p covariates, outcome
gen_data <- function(n = 100, p = 1000, model = "simple", correlated = FALSE, n_mabs = 3, strength = "strong") {
  x <- gen_x(n = n, p = p, correlated = correlated)
  ys <- gen_y(n = n, x = x, model = model, n_mabs = n_mabs, strength = strength) 
  return(list("y_bin" = ys$binary, "y_cont" = ys$continuous, "x" = x))
}