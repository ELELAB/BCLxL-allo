
path_traj=../../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate2/CHARMM22star_TIP3P/md/Mol_An/center_traj.xtc
path_tpr=../.././../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate2/CHARMM22star_TIP3P/md/md.tpr


/usr/local/plumed-2.3b/bin/plumed driver --mf_xtc ${path_traj} --plumed plumed_torsions.dat
