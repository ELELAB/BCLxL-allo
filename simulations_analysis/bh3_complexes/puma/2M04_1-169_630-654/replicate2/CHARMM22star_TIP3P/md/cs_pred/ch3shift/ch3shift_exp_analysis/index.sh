#creates index file of Bcl-xL without capping groups


path_tpr=../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate2/CHARMM22star_TIP3P/md/md.tpr



echo "a 1-2639 | a 2643\nq\n"|gmx make_ndx -f $path_tpr 
