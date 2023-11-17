# Replicating simulations

This folder contains the code necessary to replicate the simulation results in "Predicting neutralization susceptibility to combination HIV-1 monoclonal broadly neutralizing antibody regimens" by Williamson et al. (2023+). All analyses were performed using R version `4.2.2`. The simulations were run on a Windows virtual machine with 16 cores and 256GB installed RAM.

The files in this folder are:
* `00_utils.R`: generally-useful functions
* `gen_data.R`: generate a dataset, given parameters (e.g., sample size)
* `investigate_combo_timing_once.R`: run a single simulation replicate for a given set of parameters. Generate data, run prediction models, get prediction performance, return output
* `investigate_combo_timing.R`: main simulation file. Pull parameters from command line, set up parallel processing, and run the simulation a specified number of times.
* `run_all_sims.bat`: run all simulations from the Windows command line
