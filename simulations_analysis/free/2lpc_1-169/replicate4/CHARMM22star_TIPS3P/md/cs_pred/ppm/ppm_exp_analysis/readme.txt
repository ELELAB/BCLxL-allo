#This folder contains the analysis corresponding to the comparison of chemical shifts predictions from PPM
#with the experimental values (for Bcl-xL - free: 18250 (2LPC) and 18792 (2M03))

#These analysis rely on the predictions from PPM located in the folder: ../ppm_results
#Inside ../ppm_results, the results from PPM from subtrajectories with 25 ns difference are located 
#(see readme.txt inside ../ppm_results)


#The shell script ppm_exp_comparison.sh excutes the comparison of the chemical shift predictions with 
#experimental values from the corresponding chemical shift data set. For the different subtrajectories, 
#it calls the python script delta_CS.py and the output of delta_CS.py for the subtrajectories are deposited
#in the corresponding subdirectory (e.g. 25_ns_18250, 1000_ns_18792). This output comprises:

- Folder csv_summary: csv files for the output of the different atomtypes (average RMSD and reduced chi-square per atomtype) 
  and the csv with the reduced chi-sq summary (per atom) 
- Folder pdb_mapping: mapping on a pdb of the reduced-chi square (one per backbone atomtype) in blue(0)-red(3+) scale. 
  (only in folder of full trajectory)
- Folder histograms: binned distribution of the reduced chi-squared. One per backbone atom type (1 histogram for CA, 1 for CB etc.)

#The evolution of the reduced chi-sq in time for the distinct backbone atom types can be evaluated using the (crude) python script
#red_chi_sq_vs_time_plot.py. This scripts searches for the output ./cvs_summary/output_chi_squared_BB.csv in the subtrajectory
#folders, retrieves the average red chi-sq per backbone atomtype (CA,CB,N,H,HA,C) and plots the red chi-sq vs. time for all 
#backbone atom atypes. 



#-------------------------------------------#
#            delta_CS.py                    #
#-------------------------------------------#
delta_CS.py performs a comparison of predicted chemical shifts (from either PPM or ch3shift) with experimental values.
#For the particular case of PPM, the analyses and output generated have been described above. 
#delta_CS.py must be run with using python 3.7 and has the following package dependencies:
- numpy 
- sys
- re
- gromacs 
- pandas
- math
- Bio (biopandas)
- argparse
- pynmrstar 
- Matplotlib

#delta_CS.py relies on two auxiliary files contained in the folder aux_files, which must be located in the same folder as delta_CS.py:
- nomenclature.txt (nomenclature equivalencies across pdb, ch3shift and ppm)
- error_protons_ppm.txt (error for the chemical shift predictions of the different proton types by ppm)

#The python 3.7 environment on Bioinfo01 (22/04/2020) contains all these dependencies:
>> source /usr/local/envs/py37/bin/activate

#delta_CS.py assumes installation of gromacs (non mpi. You must be able to call gromacs using the alias gmx). 
#It has been tested with gromas_5.1.5, which can be source on Bioinfo01 with:
>> source /usr/local/gromacs-5.1.5/bin/GMXRC

#For help simply: 
>> python delta_CS.py -h

#An example call of delta_CS.py would look the following: 
>>  python delta_CS.py -exp BMRBENTRY -bp PATH_TO_BMDB.DAT_FROM_PPM  -ref PATH_TO_YOUR_PDB_REFERENCE -pdb_mapping yes -histograms yes

#IMPORTANT: one of the main problems when building this program was the different nomenclatures used across different force fields.
#To fix this, delta_CS.py invokes pdb2gmx and transforms the nomenclature of your reference into CHARMM27 nomenclture. 
#This implies that your reference PDB MUST NOT CONTAIN CAPPING GROUPS OR NON-CANONICAL RESIDUES NOT CONTAINED IN THIS FF. 
#Reference is used ONLY for aligning purposes, so just be sure that the numbering of your reference and your trajectory is the same.



#------------------------------------------#
#         Reproducing results              #
#------------------------------------------#

0. Source gromacs and python
>> source /usr/local/gromacs-5.1.5/bin/GMXRC
>> source /usr/local/envs/py37/bin/activate

1. Create index file to eliminate capping groups from reference:
remember that you have to remove the entire NH2 terminal capping group
ATOM   2639  NT  ARG   169      53.614  65.572  56.278  1.00  0.00           N
ATOM   2640  HT1 ARG   169      53.890  64.960  55.536  1.00  0.00           H
ATOM   2641  HT2 ARG   169      52.792  65.500  56.844  1.00  0.00           H

>> sh index.sh

2. Execute:
>> sh ppm_exp_comparison.sh

3. Plot the time evolution of the average reduced chi-sq per backbone atomtype: 
>> python red_chi_sq_vs_time_plot.py






Last update: 23/04/2020
