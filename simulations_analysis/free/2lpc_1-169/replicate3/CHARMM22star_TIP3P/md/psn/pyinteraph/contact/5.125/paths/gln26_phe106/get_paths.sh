datfilt=../../hc-graph_filt.dat
pdb=../../../../model0_cap_termodified.pdb
res1=SYSTEM-26GLN
res2=SYSTEM-106PHE

graph_analysis -r $pdb -a $datfilt -p -r1 $res1 -r2 $res2 -d -l 15
