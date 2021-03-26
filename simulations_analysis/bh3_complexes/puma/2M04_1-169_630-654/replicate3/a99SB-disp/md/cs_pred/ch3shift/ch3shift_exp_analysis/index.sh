#creates index file of Bcl-xL without capping groups
path_traj=../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate3/a99SB-disp/md/Mol_An/traj_centered.xtc
path_tpr=../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate3/a99SB-disp/md/md.tpr




echo "a 7-2645\nq\n"|gmx make_ndx -f $path_tpr 
