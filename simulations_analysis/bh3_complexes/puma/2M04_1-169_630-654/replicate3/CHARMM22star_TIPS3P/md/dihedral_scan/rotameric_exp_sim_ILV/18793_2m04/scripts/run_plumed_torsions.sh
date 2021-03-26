
path_traj=../../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate3/CHARMM22star_TIPS3P/md/Mol_An/traj_centered.xtc
path_tpr=../.././../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate3/CHARMM22star_TIPS3P/md/md.tpr


/usr/local/plumed-2.3b/bin/plumed driver --mf_xtc ${path_traj} --plumed plumed_torsions.dat
