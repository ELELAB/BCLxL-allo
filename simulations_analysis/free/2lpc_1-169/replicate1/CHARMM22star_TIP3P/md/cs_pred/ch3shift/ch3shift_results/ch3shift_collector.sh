#!/bin/bash
path_traj=../../../../../../../../../simulations/free/2lpc_1-169/replicate1/CHARMM22star_TIP3P/md/Mol_An/center_traj.xtc
path_tpr=../../../../../../../../../simulations/free/2lpc_1-169/replicate1/CHARMM22star_TIP3P/md/md.tpr

#creates index file. Removes acetylation at N terminus and NH2 at C-terminus
echo "r 1-169\nq\n"|gmx make_ndx -f ${path_tpr}


#pipeline: runs ch3dhift for the substractories 0, 25, 50, ..., 1000 ns if desired. 
for e in 1000; #0 25 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 475 500 525 500 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000; 
do
echo "19\n"|gmx trjconv -f ${path_traj} -s ${path_tpr} -b 0 -e ${e}000 -skip 100 -o traj_$e.pdb -n index.ndx #no_capping



#################################################
###########          CH3 shift          ######### ---->  ONLY FOR FULL TRAJECTORY, EXECUTION REALLY SLOW.
#################################################

#modifies command.cmd config file and remove intermediate files
sed s/TITLE/ch3_${e}/g command_template.cmd > command_1.cmd 
sed s/INPUTPDB/traj_${e}.pdb/g command_1.cmd > command_2.cmd
rm command_1.cmd
sed s/OUTPUTTXT/out_ch3shift_${e}.txt/g command_2.cmd > command.cmd
rm command_2.cmd

#executes ch3_shift
echo starting CH3shift calculations ... 
R CMD BATCH CH3Shift.R 
rm command.cmd #removes the current cmd files to leave space for the next one
rm ./traj_${e}.pdb 
cd ..
echo Done with CH3shift for t=$e ns ...


done
