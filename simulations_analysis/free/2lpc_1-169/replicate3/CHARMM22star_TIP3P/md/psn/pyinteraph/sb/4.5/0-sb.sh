#!/bin/bash

traj="../../../../../new_free/2018_analyses_for_paper/c22st_tip3p_cap/9-md/Mol_An/traj_comp_centered.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"
sb_perco=0
sb_co=4.5
ff=charmm27
sb_cg_file=charged_groups.ini

pyinteraph -v -s $gro -t $traj -r $pdb -b --sb-graph sb-graph.dat --sb-perco $sb_perco --ff-masses $ff --sb-co $sb_co --sb-cg-file $sb_cg_file

