#This folder contains the analysis corresponding to the comparison of chemical shifts predictions from CH3shift
#with the experimental values (for Bcl-xL - free: 18250 (2LPC) and 18792 (2M03))

#These analysis rely on the predictions from CH3shift located in the folder: ../ch3shift_results
#Inside ../ch3shift_results, the results from CH3shift for the full trajectory are located.
#The analysis can be run for subtrajectories as well (see readme.txt inside ../ch3shift_results for this purpose).


#The shell script ch3shift_exp_comparison.sh excutes the comparison of the chemical shift predictions with 
#experimental values from the corresponding chemical shift data sets. It calls the python script delta_CS.py
#and the output of delta_CS.py for the subtrajectories (if any) are deposited
#in the corresponding subdirectory (e.g. 25_ns_18250, 1000_ns_18792). This output comprises:
- Folder csv_summary: csv files for the output of the different atomtypes (average RMSD and reduced chi-square per atomtype) 
  and the csv file with the per-atom summary. 
- Folder pdb_mapping: mapping on a pdb of the reduced-chi square (one per atomtype, e.g. LEU CD1) in blue(0)-red(3+) scale. 
  (only in folder of full trajectory)


#-------------------------------------------#
#            Delta_CS.py                    #
#-------------------------------------------#
delta_CS.py performs a comparison of predicted chemical shifts (from either PPM or ch3shift) with experimental values.
#For the particular case of CH3shift, the analyses and output generated have been described above. 
#delta_CS.py must be run using python 3.7 and has the following package dependencies:
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

#Delta_CS.py relies on two auxiliary files contained in the folder aux_files, which must be located in the same folder as delta_CS.py:
- nomenclature.txt (nomenclature equivalencies across pdb, ch3shift and ppm)
- error_protons_ppm.txt (error for the chemical shift predictions of the different proton types by ppm)

#The python 3.7 environment on Bioinfo01 (22/04/2020) contains all these dependencies:
>> source /usr/local/envs/py37/bin/activate

#delta_CS.py assumes installation of gromacs (non mpi. You must be able to call gromacs using the alias gmx). 
#It has been tested with gromas_5.1.5, which can be source on Bioinfo01 with:
>> source /usr/local/gromacs-5.1.5/bin/GMXRC

#For help simply: 
>> python Delta_CS.py -h

#An example call of Delta_CS.py would look the following: 
>>  python delta_CS.py -exp BMRBENTRY -ch3 PATH_TO_BMDB.DAT_FROM_CH3shift -ref PATH_TO_YOUR_PDB_REFERENCE -pdb_mapping yes -histograms no

#IMPORTANT: one of the main problems when building this program was the different nomenclatures used across different force fields.
#To fix this, delta_CS.py invokes pdb2gmx and transforms the nomenclature of your reference into CHARMM27 nomenclture. 
#This implies that your reference PDB MUST NOT CONTAIN CAPPING GROUPS OR NON-CANONICAL RESIDUES NOT CONTAINED IN THIS FF. 
#Reference is used ONLY for aligning purposes, so just be sure that the numbering of your reference and your trajectory is the same.



#------------------------------------------#
#         Reproducing results              #
#------------------------------------------#

0. Source gromacs and source the python environment.
>> source /usr/local/gromacs-5.1.5/bin/GMXRC
>> source /usur/local/envs/py37/bin/activate

1. Create index file to eliminate capping groups from reference:
>> sh index.sh

2. Execute:
>> sh ch3shift_exp_comparison.sh








