traj="../../../../new_complexes_capped/bclxl-puma/9-md/Mol_An/traj_centered.xtc"
tpr="../../../../new_complexes_capped/bclxl-puma/9-md/md.tpr"
ndx="A-E153_A-R164_sidechains.ndx"

g_dist -f $traj -s $tpr -n $ndx -o dist_E153_R164_sidechains.xvg <<eof
0
1
eof
