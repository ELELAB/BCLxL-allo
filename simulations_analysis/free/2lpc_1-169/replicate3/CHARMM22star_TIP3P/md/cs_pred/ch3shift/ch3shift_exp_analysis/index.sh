#creates index file of Bcl-xL without capping groups (NH2 at C terminus, i.e. deleting 2640, 2641, and 2642)


echo "a 1-2639 | a 2643 \nq\n"|gmx make_ndx -f ../../../../../../../../../simulations/free/2lpc_1-169/replicate3/CHARMM22star_TIP3P/md/md.tpr
