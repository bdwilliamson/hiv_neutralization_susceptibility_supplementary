# simulation to check if combining bnAb susceptibility before or after making predictions is better
# vary with sample size

# load required functions and packages -----------------------------------------
library("dplyr")
library("tidyr")
library("data.table")
library("cvAUC")
library("glmnet")
library("parallel")
library("rprojroot")
library("optparse")

this_path <- normalizePath(".", mustWork = FALSE)
proj_root <- rprojroot::find_root_file(criterion = ".projectile", path = this_path)

source(paste0(proj_root, "/code/sims/gen_data.R"))
source(paste0(proj_root, "/code/sims/investigate_combo_timing_once.R"))
source(paste0(proj_root, "/code/R/00_utils.R"))

# get command-line args --------------------------------------------------------
parser <- OptionParser()
parser <- add_option(parser, "--n", type = "integer", default = 500,
                     help = "sample size")
parser <- add_option(parser, "--strength", default = "strong", 
                     help = "strength of association between covariates and outcome")
parser <- add_option(parser, "--nreps-total", type = "integer", default = 2500,
                     help = "Number of total replicates")
parser <- add_option(parser, "--nreps-per-job", type = "integer", default = 1,
                     help = "Number of reps to run per job")
args <- parse_args(parser, convert_hyphens_to_underscores = TRUE)
print(args)

output_dir <- paste0(proj_root, "/results/sims/")
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

all_analyses <- expand.grid(
  n = seq(100, 1000, 100)
)
n_jobs <- args$nreps_total / args$nreps_per_job
job_id <- which((args$n == all_analyses$n))

# set up arguments for data generation -----------------------------------------
model <- "simple"
p <- 1000
n_mabs <- 3

# set up parallelization -------------------------------------------------------
num_cores <- parallel::detectCores()
cl <- parallel::makePSOCKcluster(num_cores)
parallel::clusterEvalQ(cl, library("glmnet"))
parallel::clusterExport(cl, varlist = c(
  "gen_x", "gen_y", "gen_data", "investigate_combo_timing_once", "cvR2",
  "additive_combo_method", "get_cv_preds", "model", "p", "n_mabs", "args",
  "n_jobs"
))
seed <- job_id * 1000
set.seed(seed)
clusterSetRNGStream(cl = cl, iseed = seed)
# run the simulation
cat("Running analysis\n")
start <- Sys.time()
output_list <- parallel::parLapply(
  cl = cl, X = as.list(1:args$nreps_per_job),
  fun = function(i) {
    investigate_combo_timing_once(
      mc_id = i + (n_jobs - 1) * args$nreps_per_job, n = args$n, p = p, 
      model = model, n_mabs = n_mabs, correlated = FALSE, strength = args$strength
    )}
)
end <- Sys.time()
cat("Elapsed time: ", format(end - start), "\n")
output <- data.table::rbindlist(output_list)
parallel::stopCluster(cl)

saveRDS(output, file = paste0(output_dir, "/output_n", args$n, "_", args$strength,
                              ".rds"))