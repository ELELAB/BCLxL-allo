#!/bin/bash

. /usr/local/gromacs-5.1.5/bin/GMXRC

traj="../../../new_complexes_capped_small_box/bclxl-puma/9-md/Mol_An/traj_centered.xtc"
tpr="../../../new_complexes_capped_small_box/bclxl-puma/9-md/md.tpr"
out="model0.pdb"

echo "1\n" | gmx trjconv -f $traj -s $tpr -o $out -b 0 -e 1
