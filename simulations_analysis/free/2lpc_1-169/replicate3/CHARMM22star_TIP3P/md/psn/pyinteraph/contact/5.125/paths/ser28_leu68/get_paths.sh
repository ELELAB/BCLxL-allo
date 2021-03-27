datfilt=../../hc-graph_filt.dat
pdb=../../../../model0_cap_termodified.pdb
res1=SYSTEM-28SER
res2=SYSTEM-68LEU

graph_analysis -r $pdb -a $datfilt -p -r1 $res1 -r2 $res2 -d -l 15
