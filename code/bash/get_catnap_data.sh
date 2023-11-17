#!/bin/bash 

# all of the individual bnAbs from the combinations from Williamson et al. (2023+) (may need more)
# plus nAbs from the SLAPNAP software paper: all combos with > 125 evaluable combination sequences and individual sequences
ALL_NABS=("VRC07-523-LS" "PGDM1400" "VRC26.25" "10-1074" "PGT121" \
          "10E8" "PG9" "VRC07" "3BNC117" "PGT128" "10-1074" "NC37" "BG1" "BG18" \
          "10E8+PG9+VRC07" "10E8+3BNC117+PG9" \
          "10E8+PG9+PGT128" "10-1074+10E8" "10-1074+10E8+3BNC117" \
          "10-1074+3BNC117" \
          "10-1074+3BNC117+PG9" "10-1074+PG9" "10E8+3BNC117" \
          "10E8+PGT128+VRC07" "3BNC117+PG9" \
          "BG1+BG18+NC37" "PG9+PGT128" \
          "PG9+PGT128+VRC07" "PG9+VRC07" "PGT128+VRC07")
# ALL_NABS=("VRC07-523-LS+PGDM1400" \
#           "VRC07-523-LS+VRC26.25" \
#           "VRC07-523-LS+10-1074" \
#           "VRC07-523-LS+PGT121" \
#           "VRC07-523-LS+PGT121+PGDM1400")
# combo_methods=("additive" "Bliss-Hill")

# loop through and grab catnap data; note it will get all outcomes
# for combo_method in ${combo_methods[@]}; 
# do
#   mkdir -p ./data/${combo_method}
# SLAPNAP_NAB=$(echo $NAB | sed 's/+/;/g' | sed 's/ //g')
for NAB in ${ALL_NABS[@]}
do
  DATA_NAB=$(echo $NAB | sed 's/ //g')

  sudo docker run \
    -v "$(pwd)"/data/:/home/output \
    -e nab=$DATA_NAB \
    -e outcomes="ic50;ic80" \
    -e return="data" \
    slapnap/slapnap
done
# done