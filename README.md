# Predicting neutralization susceptibility to combination HIV-1 bnAb regimens

This repository contains the code necessary to replicate the results in "Predicting neutralization susceptibility to combination HIV-1 monoclonal broadly neutralizing antibody regimens" by Williamson et al. (2023+). All analyses were performed using R version >= `4.2.2`; the data analyses used R version `4.3.1` while the simulations used R version `4.2.2`.

The `code` folder has three further subfolders:
1. `bash`: contains shell scripts (`*.sh`) and Windows batch files (`*.bat`) that run certain R scripts
2. `R`: the R code to run analyses
3. `sims`: code to replicate simulations

## bash
* `bash/get_catnap_data.sh` obtains all of the `data/*.csv` files from CATNAP
* `bash/get_all_neut_data.bat`: wrangles the CATNAP data for all individual bnAbs and bnAb regimens into the format we need for the analyses
* `bash/get_all_preds.bat`: gets all predictions for all bnAb regimens, based on all outcomes and combo methods/times

## R
* `R/00_utils.R`: contains generally-useful functions
* `R/01_combo_neut.R`: wrangle the data into the format we need for the analysis. Called by `bash/get_all_neut_data.bat`
* `R/02_get_predictions.R`: gets all predictions for a given bnAb regimen. Called by `bash/get_all_preds.bat`
* `R/03_prediction_performance.R`: obtains prediction performance for a given bnAb regimen.
* `R/04_analysis_plots.R`: plots of results and data for a given bnAb regimen.

## sims

Contains code to replicate the simulations.