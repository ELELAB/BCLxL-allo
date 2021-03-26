PATHdata='../../../../../../../../simulations/free/2lpc_1-169/replicate4/a99SB-disp/md/Mol_An2'
xtc=center_traj
/usr/local/plumed-2.3b/bin/plumed driver --mf_xtc $PATHdata/$xtc.xtc --plumed plumed_torsions.dat
