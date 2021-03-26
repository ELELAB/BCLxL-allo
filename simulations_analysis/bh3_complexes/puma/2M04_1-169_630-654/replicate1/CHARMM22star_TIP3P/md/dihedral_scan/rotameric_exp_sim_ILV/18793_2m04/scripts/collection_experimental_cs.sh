source /usr/local/envs/py37/bin/activate



#creates reference
path_traj=../../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate1/CHARMM22star_TIP3P/md/Mol_An/traj_centered.xtc
path_tpr=../../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate1/CHARMM22star_TIP3P/md/md.tpr

echo "a 1-2639\nq\n"| gmx make_ndx -f ${path_tpr}
echo "19\n"|gmx trjconv -f ${path_traj} -s ${path_tpr} -b 0 -e 0 -o reference.pdb -n index.ndx




bmrbentry=18793
#cp ../delta_CS.py .
#cp -r ../aux_files .
python delta_CS.py -exp $bmrbentry  -ref reference.pdb 

rm ./#*
