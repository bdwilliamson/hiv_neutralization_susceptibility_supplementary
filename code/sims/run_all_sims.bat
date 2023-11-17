:: run all scenarios for all non-TMLE estimators
@echo off
setlocal enabledelayedexpansion

:: set R version for your machine
:: set r_version=4.3.1
:: set r_version=4.2.1
For /F "Delims=" %%0 In ('where /r "C:\Program Files\R" Rscript.exe') do set scriptpath="%%~0"
echo Using R executable !scriptpath!

:: set up arguments to pass
set nreps=2500
:: set strength=weak
:: send output to directory
set outdir=%cd%\rout

if not exist %outdir% mkdir %outdir%

:: run the simulation
:: outer loop over strength
:: inner loop over sample size
for %%S in ("strong" "weak") do (
  for %%N in (100 200 300 400 500 600 700 800 900 1000) do (
    set this_outfile=!outdir!\output_n%%N_%%S.out
    echo Running n = %%N, %%S
    !scriptpath! investigate_combo_timing.R --nreps-total !nreps! --n %%N --strength %%S --nreps-per-job !nreps! 1>!this_outfile! 2>&1
  )
)

