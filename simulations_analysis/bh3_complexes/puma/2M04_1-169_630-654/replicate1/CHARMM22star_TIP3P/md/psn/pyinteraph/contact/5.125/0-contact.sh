#!/bin/bash
traj="../../../new_complexes_large_box/complexes/puma_3.5nm/9-md/traj_centered.xtc"
gro="../model0.pdb"
pdb="../model0.pdb"


pyinteraph -s $gro -t $traj -r $pdb -f --hc-graph hc-graph.dat --ff-masses charmm27 --hc-co 5.125 --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,ARG,HIS,LYS,GLU,ASP,ASN,GLN,SER,THR,CYS

