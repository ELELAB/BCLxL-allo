#!/bin/bash
traj="../../../../../new_free/2018_analyses_for_paper/c22st_tip3p_cap/9-md/Mol_An/traj_comp_centered.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"
hc_co=5.125
hc_perco=0
ff=charmm27

pyinteraph -s $gro -t $traj -r $pdb -f --hc-graph hc-graph.dat --ff-masses $ff --hc-co $hc_co --hc-perco $hc_perco --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,ARG,HIS,LYS,GLU,ASP,ASN,GLN,SER,THR,CYS

