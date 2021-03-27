traj="../../../../new_complexes_large_box/complexes/puma_3.5nm/9-md/traj_centered.xtc"
tpr="../../../../new_complexes_large_box/complexes/puma_3.5nm/9-md/md.tpr"
ndx="A-H73_A-T75_sidechains.ndx"

g_dist -f $traj -s $tpr -n $ndx -o dist_H73_T75_sidechains.xvg <<eof
0
1
eof
