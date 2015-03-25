This directory contains the input data, protocols and output model for the
modeling of the PhoQ homodimer, using cysteine crosslinking and
multi-state Bayesian modeling in [IMP](http://integrativemodeling.org).

## Directories:

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

To run a modeling script, just change into its directory and run it from the
command line, e.g.
 - `cd modeling/2Y20_MC_2_FULL`
 - `./HK_model_ENSEMBLE_REM.py` (on a single processor; prepend `mpirun -np 4` or similar if you built IMP with MPI support)

Two output files will be produced: a log file containing statistics, and a
trajectory of model structures in [RMF](http://integrativemodeling.org/rmf/)
format. The latter can be viewed in [UCSF Chimera](http://www.cgl.ucsf.edu/chimera/).

## Information

_Author(s)_: Max Bonomi, Riccardo Pellarin, Ben Webb

_Date_: March 18th, 2015

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Last known good IMP version_: [![build info](https://salilab.org/imp/systems/?sysstat=12)](http://salilab.org/imp/systems/)

_Testable_: Yes.

_Parallelizeable_: Yes

_Publications_:
 - K. Molnar, M. Bonomi, R. Pellarin, G. Clinthorne, G. Gonzalez, S. Goldberg, M. Goulian, A. Sali, W. DeGrado. [Cys-Scanning Disulfide Crosslinking and Bayesian Modeling Probe the Transmembrane Signaling Mechanism of the Histidine Kinase, PhoQ](http://www.ncbi.nlm.nih.gov/pubmed/25087511), Structure 22, 1239-1251, 2014.
