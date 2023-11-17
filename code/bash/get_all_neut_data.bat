:: get all combo neutralization data
@echo off
setlocal enabledelayedexpansion

set outfile=%cd%\rout
if not exist %outfile% mkdir %outfile%

set r_version=4.3.1

for %%B in ("VRC07-523-LS+PGDM1400" "VRC07-523-LS+VRC26.25" "VRC07-523-LS+10-1074" "VRC07-523-LS+PGT121" "VRC07-523-LS+PGT121+PGDM1400" "10E8+PG9+VRC07" "10E8+3BNC117+PG9" "10E8+PG9+PGT128" "10-1074+10E8" "10-1074+10E8+3BNC117" "10-1074+3BNC117" "10-1074+3BNC117+PG9" "10-1074+PG9" "10E8+3BNC117" "10E8+PGT128+VRC07" "3BNC117+PG9" "BG1+BG18+NC37" "PG9+PGT128" "PG9+PGT128+VRC07" "PG9+VRC07" "PGT128+VRC07") do (
  set this_outfile=%outfile%\output_%%B.out 
  echo outfile: !this_outfile!
  "C:\Program Files\R\R-!r_version!\bin\Rscript.exe" R\01_combo_neut.R --bnab %%B 1>!this_outfile! 2>&1
)
