datfilt=../../hc-graph_filt.dat
pdb=../../../../model0.pdb
res1=SYSTEM-23SER
res2=SYSTEM-73HIS

graph_analysis -r $pdb -a $datfilt -p -r1 $res1 -r2 $res2 -d -l 15
