PATHdata='../../../../../../../../simulations/free/2lpc_1-169/replicate3/CHARMM22star_TIP3P/md/Mol_An'
xtc=center_traj.xtc
/usr/local/plumed-2.3b/bin/plumed driver --mf_xtc $PATHdata/$xtc --plumed plumed_torsions.dat
