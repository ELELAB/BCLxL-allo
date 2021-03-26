#creates index file of Bcl-xL without capping groups (NH2 at C terminus and acetylation at N terminus, i.e. ri 2-170 [or r 1-169])


echo "ri 2-170\nq\n"|gmx make_ndx -f ../../../../../../../../../../../simulations/bh3_bcl/bclxl/free/2lpc_1-169/replicate4/a99SB-disp/md/md.tpr
