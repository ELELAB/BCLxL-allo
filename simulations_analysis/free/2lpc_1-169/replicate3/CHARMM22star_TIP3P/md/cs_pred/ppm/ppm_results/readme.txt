#This folder contains the results by PPM, i.e. chemical shift predictions for backbone atoms and siide-chain hydrogens.
#PPM takes as input a centered PDB trajectory (alghouth PPM compensates for potential PBC defects).

#The script ppm_collector.sh assumes the presence of the executable ppm_stable in the folder
#The paths of the TPR and XTC files (from which PDBs are created) are defined at the beginning of the script ppm_collector.sh:
>> path_traj=/../../../../../../../../../simulations/free/2lpc_1-169/replicate4/a99SB-disp/md/Mol_An2/center_traj.xtc
>> path_tpr=../../../../../../../../../simulations/free/2lpc_1-169/replicate4/a99SB-disp/md/md.tpr

#The script ppm_collector.sh collects the chemical shift predictions for subtrajectories of 25 ns lenght difference, 
#i.e. for 0ns, 25ns, 50ns, ..., 1000ns. The script creates the intermediate subtrajectories and deletes them as it progresses
#in the loop. 

#The input trajectory was processed by removing the capping groups. The line to create the index file is at the beginning of the ppm_collector.sh
#script. 

#To reproduce results, simply run: 
>> sh ppm_collector.sh
