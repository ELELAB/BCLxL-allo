#!/bin/bash
traj="../../../new_free/2018_analyses_for_paper/c22st_tips3p/9-md_replica1_box2nm/Mol_An/traj_comp_centered.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"


pyinteraph -s $gro -t $traj -r $pdb -f --hc-graph hc-graph.dat --ff-masses charmm27 --hc-co 5.125 --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,ARG,HIS,LYS,GLU,ASP,ASN,GLN,SER,THR,CYS

