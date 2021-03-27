#!/bin/bash

dat="hc-graph.dat"
datfilt="hc-graph_filt.dat"
gro="../../../../../new_free/2018_analyses_for_paper/c22st_tip3p/9-md_replica2_box2nm/Mol_An/model0.pdb"
pdb="../../../../../new_free/2018_analyses_for_paper/c22st_tip3p/9-md_replica2_box2nm/Mol_An/model0.pdb"

# Filter the graph
filter_graph -d $dat -o $datfilt -t 20.0

# Write the hubs
graph_analysis -a $datfilt -r $pdb -u -ub hubs.pdb -k 3

# Write the connected components
graph_analysis -a $datfilt -r $pdb -c -cb con_comp.pdb

