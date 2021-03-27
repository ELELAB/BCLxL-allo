#!/bin/bash

traj="../../../../../new_complexes_capped_small_box/bclxl-puma/9-md/Mol_An/traj_centered.xtc"
gro="../../model0.pdb"
pdb="../../model0.pdb"
hc_perco=0
hc_co=5.125
ff=charmm27

pyinteraph -s $gro -t $traj -r $pdb -f --hc-graph hc-graph.dat --hc-perco $hc_perco --ff-masses $ff --hc-co $hc_co --hc-residues ALA,ILE,LEU,VAL,PHE,TRP,TYR,MET,PRO,ARG,HIS,LYS,GLU,ASP,ASN,GLN,SER,THR,CYS

