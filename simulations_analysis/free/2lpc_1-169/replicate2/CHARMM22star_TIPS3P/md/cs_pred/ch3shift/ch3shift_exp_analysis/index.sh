#creates index file of Bcl-xL without capping groups (NH2 at C terminus and acetylation at N terminus, i.e. ri 2-170 [or r 1-169])


echo "r 1-169\nq\n"|gmx make_ndx -f ../../../../../../../../../simulations/free/2lpc_1-169/replicate2/CHARMM22star_TIPS3P/md/md.tpr
