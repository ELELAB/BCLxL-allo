#!/bin/bash

dat="hc-graph.dat"
datfilt="hc-graph_filt.dat"
gro="../../model0.pdb"
pdb="../../model0.pdb"

# Filter the graph
filter_graph -d $dat -o $datfilt -t 20.0

# Write the hubs
graph_analysis -a $datfilt -r $pdb -u -ub hubs.pdb -k 3

# Write the connected components
graph_analysis -a $datfilt -r $pdb -c -cb con_comp.pdb

