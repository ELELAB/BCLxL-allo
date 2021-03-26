# BCLxL-allo
Computational Biology Laboratory, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark
Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark


Repository associated to the publication:

Bcl-xL dynamics under the lens of protein structure networks. Valentina Sora, Dionisio Sanchez, Elena Papaleo*,J Phys Chem B, under revision
a early pre-print version in biorxiv: doi: https://doi.org/10.1101/574699

contacts for repository:elenap-at-cancer.dk

The repository contains the input data for modeling and simulations, outputs from analyses and associated scripts to reproduce our data. The the MD trajectories, mutational scans have been deposited in OSF due to space limitation. They can be downloaded from here:

https://osf.io/rnc6d/

In simulations_analysis folders with filtered tpr or xtc files do not contain the outputs due to space limitation. Everything is reproducible through the scripts.

A more detailed README file is in each subfolder

pdbs/ -> preprocessing and cleaning steps for the PDB files
simulations_analysis -> all the analyses done on free and PUMA-bound simulations
mutatex -> foldx scans a saturation for free protein or complexes


In the paper we also made usage of the code deposited in: https://github.com/ELELAB/delta_cs, to compare experimental and predicted chemical shifts. 

Please cite our publication if you use the material in this repository

Software versions: gromacs , foldx-suite



