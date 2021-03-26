
#!/bin/bash
path_traj=../../../../../../../../../simulations/free/2lpc_1-169/replicate3/CHARMM22star_TIP3P/md/Mol_An/center_traj.xtc
path_tpr=../../../../../../../../../simulations/free/2lpc_1-169/replicate3/CHARMM22star_TIP3P/md/md.tpr

#creates index file.
echo "r 1-169\nq\n"|gmx  make_ndx -f ${path_tpr}


#pipeline: runs ppm for the substractories 0, 25, 50, ..., 1000 ns and ch3shift only for the complete trajectory (1000 ns)
for e in 0 25 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 475 500 525 550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000; 
do
echo "19\n"|gmx  trjconv -f ${path_traj} -s ${path_tpr} -b 0 -e ${e}000 -skip 100 -o traj_$e.pdb -n index.ndx #no_capping


###################################################
############            PPM           ############
##################################################


mkdir ppm_${e}
cd ppm_${e}
echo Starting with ppm calculations ... 
../ppm_stable -pdb ../traj_${e}.pdb #> ppm_${e}_output.log &

 # However, it does not work when the sequences are not identical. It did not work with our trajs. 

echo Done with PPM for t=$e ns ...
cd ..

rm ./#*
#remove trajectory
echo Removing created trajectory ... 
rm traj_${e}.pdb

done



