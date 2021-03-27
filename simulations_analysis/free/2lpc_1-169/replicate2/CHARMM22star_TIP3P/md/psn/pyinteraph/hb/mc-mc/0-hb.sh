#!/bin/bash

traj="../../../../../new_free/2018_analyses_for_paper/c22st_tip3p/9-md_replica2_box2nm/Mol_An/traj_comp_centered.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"
hb_perco=0
ff="charmm27"
hb_class="mc-mc"
hb_ad_file="hydrogen_bonds.ini"

pyinteraph -v -s $gro -t $traj -r $pdb -y --hb-graph hb-graph.dat --hb-perco $hb_perco --ff-masses $ff --hb-class $hb_class --hb-ad-file $hb_ad_file

~                                                                               
