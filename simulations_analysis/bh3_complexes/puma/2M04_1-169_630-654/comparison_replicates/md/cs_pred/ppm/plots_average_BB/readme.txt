This folder contains the plots comparing the results of the average reduced chi squared for backbone atoms computed from 
the trajectories of the different replicates with the goal of assessing the degree of agreement of each model with 
the available chemical shift data.

For the trajectories with PUMA, 2 different visualizations are made available: 
1. Average reduced chi squared for the backbone atoms. 
		--> plot_comparison.py
2. Average reduced chi squared for the backbone atoms, excluding 5 outliers for CA and CB driving the metric to extreme values, 
   potentially not accounting for the real level of agreement. The motivation to do so is that only those 5 atoms accounted for 
   more than half of the total mass of the reduced chi squared metric. Moreover, for those "conflictive residues", the bad 
   chi-squared are not present in all backbone atom types (i.e. when an anomally high average reduced chi-squared was present
   for a CA, for instance, this trend was not observed in the rest of the bb atom types: HA, H, CB, etc.) 
		--> plot_comparison.py (creates filtered statistics under "filtered_csv")
		The residues filtered out are below. In brackets, the residues numbering in Uniprot (paper)-
		Filtered CB: F103 (143), H73 (113), F65 (105), T75 (115), E118 (158)
		Filtered CA: D55 (95), F104 (144), L59 (99), W129 (169),  M119 (159)  


###############################
##       REPRODUCE        #####
###############################

0. source /usr/local/envs/py37/bin/activate
1. python plot_comparison.py 
2. python plot_comparison_filtered.py
