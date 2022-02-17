# BCLxL-allo
Computational Biology Laboratory, Danish Cancer Society Research Center, 2100, Copenhagen, Denmark
Cancer Systems Biology, Health and Technology Department, Section for Bioinformatics, 2800, Lyngby, Denmark


Repository associated to the publication:

Bcl-xL dynamics under the lens of protein structure networks. Valentina Sora, Dionisio Sanchez, Elena Papaleo*, J Phys Chem B
. 2021 May 6;125(17):4308-4320. doi: 10.1021/acs.jpcb.0c11562. Epub 2021 Apr 13.

contacts for repository:elenap-at-cancer.dk

The repository contains the input data for modeling and simulations, outputs from analyses and associated scripts to reproduce our data. The the MD trajectories have been deposited in OSF due to space limitation. They can be downloaded from here:

https://osf.io/rnc6d/

In simulations_analysis folders with filtered tpr or xtc files do not contain the outputs due to space limitation. Everything is reproducible through the scripts.

A more detailed README file is in each subfolder

simulations_analysis -> all the analyses done on free and PUMA-bound simulations

mutatex -> foldx scans a saturation for free protein or complexes


In the paper we also made usage of the code deposited in: 

https://github.com/ELELAB/delta_cs to compare experimental and predicted chemical shifts. 

https://github.com/ELELAB/psntools to carry out centrality measurements and other additional PSN analyses
which complement the ones available in PyInteraph 


Please cite our publication if you use the material in this repository

Software versions: gromacs version 4 and 5, foldx-suite4, pyinteraph

N.B. The replicates folders are numbered from 1 to 4 but for some of them two different force fields or two different solvent models have been used for a total count of 7 replicates (i.e. independent trajectories) of the free BCL-xL and four of the complex BCL-xL:PUMA



