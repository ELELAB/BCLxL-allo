#!/bin/bash

. /usr/local/gromacs-5.1.5/bin/GMXRC

traj="../../../new_free/2018_analyses_for_paper/c22st_tip3p/9-md_replica2_box2nm/Mol_An/traj_comp_centered.xtc"
tpr="../../../new_free/2018_analyses_for_paper/c22st_tip3p/9-md_replica2_box2nm/md.tpr"
out="model0.pdb"

echo "1\n" | gmx trjconv -f $traj -s $tpr -o $out -b 0 -e 1
