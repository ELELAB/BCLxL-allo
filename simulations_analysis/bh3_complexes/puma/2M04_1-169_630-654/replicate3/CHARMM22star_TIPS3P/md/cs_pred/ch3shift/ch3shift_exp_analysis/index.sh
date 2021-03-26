#creates index file of Bcl-xL without capping groups


path_tpr=../../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate3/CHARMM22star_TIPS3P/md/md.tpr



echo "a 7-2644 | a 2648\nq\n"|gmx make_ndx -f $path_tpr 
