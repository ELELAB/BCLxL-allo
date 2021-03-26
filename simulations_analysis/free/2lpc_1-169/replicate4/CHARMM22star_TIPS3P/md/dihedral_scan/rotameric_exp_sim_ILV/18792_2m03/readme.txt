#This folder contains the analysis corresponding to the comparison of the rotameric side-chain 
#populations of Ile, Val and Leu from chemical shifts and simulations. This analysis provides further insight into
#agreement of side-chain dynamics of MD simulations with experimental data. 

#The rotameric populations of Val, Ile and Leu side-chains are defined according to the chi1 (Val), chi2(Ile) and chi2(Leu)
#The atoms defining these dihedrals can be consulted in http://www.ccp14.ac.uk/ccp/web-mirrors/garlic/garlic/commands/dihedrals.html
#These rotameric states are: 
 - Trans: from -180, to -120 and from 120 to 180
 - Gauche +: from 0 to 120 
 - Gauche -: from -120 to 0


#The side-chain rotameric populations for Val, Ile and Leu can be estimated from chemical shifts, following the following methods 
# - Leu: from the difference of the chemical shifts of the atoms CD1 and CD2, trans and g+ populations (of chi2) can be estimated (Mulder 2009) 
	Pt = [(δCD1-δCD2)+ 5 ppm]/10 ppm, Pg+ = 1-Pt (only g+ and t+ populated)
# - Val: from the chemical shift pair of CG1 and CG2, the trans, g+ and g- populations (of chi1) can be estimated (Hansen and Kay 2011).
	More complex. See reference.
# - Ile: from the chemical shift CD1, the trans and g- populations (of chi1) can be estimated (Hansen 2010)
	Pg ={
	1 for δ(CD1) < 9.3ppm
	0 for δ(CD1) > 14.8ppm
	else
	[14.8ppm - δ(CD1)]/5.5ppm
#These results are compiled in the web-server application SideR by Hansen's group: https://www.ucl.ac.uk/hansen-lab/sider_dfh.html


#--------------------------------------------------------------------------------------------------------#
#                           FOLDER CONTENTS                                                              #
#--------------------------------------------------------------------------------------------------------#
- A directory called "scripts" with all the shell and python scripts involved in the pipeline, including input and intermediate 
  output files. These scripts are described below under "SCRIPTS AND DEPENDENCIES"
- A .csv summary called "table_comparison_exp_sim_18792.csv", table compiling the rotameric populations (experimental and sim)
  for all Ile, Leu and Val residues in our structure (Bcl-xL)
- A folder called "plots" with the plots comparing simulation and experiment for Val, Ile and Leu residues.  



#--------------------------------------------------------------------------------------------------------#
#                           SCRIPTS AND DEPENDENCIES                                                     #
#--------------------------------------------------------------------------------------------------------#
#All scripts as well as their input files are located within the folder scripts. Intermediate outputs of the pipeline are stored
#that folder too.

-  delta_CS.py: used to retrieve the chemical shifts from the BMRB entry (see:  ../../../cs_pred/ppm/ppm_exp_analysis/
	for more information on dependencies)
-  collection_experimental_cs.sh: shell script with the command line to execute delta_CS.py for retrieval of chemical shifts. 
-  collection_input_sideR.py: python script for the creation of the input for sideR (python 3, pandas, os, NumPy)
- plotter_torsions.py: python script to merge results from sideR and the simulations and to generate the output plots (python 3,
	NumPy, Pandas, Matplotlib, Re, Biopandas, os).  
-  create_plumed_torsions.py: a python script that creates (from the plumed template plumed_template.dat) a plumed file (plumed_torsions.dat)
	with the commands to compute and dump to a text file (torsion all) the torsion angles of interest. 




#--------------------------------------------------------------------------------------------------------#
#                           REPRODUCTION                                                                 #
#--------------------------------------------------------------------------------------------------------#


Execute the analysis described below within the folder "scripts" (cd scripts). To reproduce the analysis: 

0. source Python and GROMACS
>> source /usr/local/envs/py37/bin/activate
>> source /usr/local/gromacs-5.1.5/bin/GMXRC

1. Collect the chemical shifts from the corresponding BMRB  entry and dump to a csv file called experimental.csv. 
   This step generates the file experimental.csv. ATTENTION: check the field bmrbentry in the .sh file.
   This step is done with the shell script "collection_experimental_cs.sh", that calls "delta_cs.py" .
>> sh collection_experimental_cs.sh 

2. Run the python script "collection_input_sideR.py" to retrieve the relevant chemical shifts entries and produce an input .txt
   file compatible with sideR. This step generates the file input_file_sider.txt 
>> python  collection_input_sideR.py

3. Copy the content of the file input_file_sider.txt (ATTENTION: be careful that you don't add a blank line at the end or sideR will complain)
   and insert in the input window https://www.ucl.ac.uk/hansen-lab/sider_dfh.html. Submit and wait for the output. 

4. Copy THE ENTIRE page prompted as output and paste it into a file called output_sider.txt (within the scripts folder)


5. (Done already) Create a file called "plumed_template.dat" (only with the entries MOLINFO and WHOLEMOLECULES). E.g.:

'''			#No need of restart since -append in gromacs automatically sets same behaviour as RESTART
			MOLINFO STRUCTURE=reference.pdb MOLTYPE=protein

			#Wholemolecules entry (BCL-XL)
			WHOLEMOLECULES ENTITY0=1-2648 # For fixing PBC problems                                       '''


6. Create the file "plumed_torsions.dat" to compute the chi angles of interest of Val, Leu and Ile residues from our trajectory. 
>> python create_plumed_torsions.py


7. Run plumed on the trajectory using the file torsion_all. This is done in the shell script "run_plumed_torsions.sh"
>> sh run_plumed_torsions.sh

8. Run the python script "plotter_torsions.py", to compare the dihedrals populations. ATTENTION: change the field dataset of the python script if necessary. 
>> python plotter_torsion.py 





LAST UPDATED: 25/04/2020
