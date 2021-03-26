#creates index file of Bcl-xL without capping groups


path_tpr=../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate1/CHARMM22star_TIP3P/md/md.tpr



echo "r 1-169\nq\n"|gmx make_ndx -f $path_tpr 
