datfilt=../../hc-graph_filt.dat
pdb=../../../../model0.pdb
res1=A-22TYR
res2=A-68LEU

graph_analysis -r $pdb -a $datfilt -p -r1 $res1 -r2 $res2 -d -l 20
