source /usr/local/envs/py37/bin/activate



#creates reference
path_traj=../../../../../../../../../../../../simulations/bh3_bcl/bclxl/free/2lpc_1-169/replicate4/a99SB-disp/md/Mol_An/center_traj.xtc
path_tpr=../../../../../../../../../../../../simulations/bh3_bcl/bclxl/free/2lpc_1-169/replicate4/a99SB-disp/md/md.tpr


echo "1\n"|gmx trjconv -f ${path_traj} -s ${path_tpr} -b 0 -e 0 -o reference.pdb




bmrbentry=18792
cp ../delta_CS.py .
cp -r ../aux_files .
python delta_CS.py -exp $bmrbentry  -ref reference.pdb 

rm ./#*
