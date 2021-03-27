traj="../../../../new_free/2018_analyses_for_paper/c22st_tips3p/9-md_replica1_box2nm/Mol_An/traj_comp_centered.xtc"
tpr="../../../../new_free/2018_analyses_for_paper/c22st_tips3p/9-md_replica1_box2nm/md.tpr"
ndx="H73_T75_sidechains.ndx "

g_dist -f $traj -s $tpr -n $ndx -o dist_H73_T75_sidechains.xvg <<eof
0
1
eof
