This directory contains the models for the structure of the phoq homodimer,
using cysteine crosslinking and 2-state Bayesian modeling.

`solution_ensemble`:
all the 2-state models within the same cluster of similar structures.
Each model is made of two structures, labeled by number 1 and 2.

`localization_densities`:
all the calculated localization densities, for given domains, and for the
two states.

`ucsf_chimera_files`:
files used to generate the figures in the paper.

`modeling`:
scripts for modeling the system:

 - `2Y20_MC_1_FULL` 1-state modeling
 - `2Y20_MC_2_FULL` 2-state modeling
 - `2Y20_MC_2_FULL_90_PERCENT_DATA_1` 2-state modeling with 90% data jackknifing
 - `2Y20_MC_2_FULL_90_PERCENT_DATA_2` 2-state modeling with 90% data jackknifing
 - `2Y20_MC_2_FULL_95_PERCENT_DATA_1` 2-state modeling with 95% data jackknifing
 - `2Y20_MC_2_FULL_95_PERCENT_DATA_2` 2-state modeling with 95% data jackknifing
 - `2Y20_MC_3_FULL` 3-state modeling
