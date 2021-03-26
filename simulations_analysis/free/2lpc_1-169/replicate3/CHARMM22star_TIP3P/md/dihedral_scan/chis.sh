source /usr/local/gromacs-2019.4/bin/GMXRC  #do not use gromacs5 for this due to a bug which gives segmentation fault.
PATHdata='../../../../../../../../simulations/free/2lpc_1-169/replicate3/CHARMM22star_TIP3P/md/Mol_An'
gro=update
trj=center_traj
mkdir chis
cd chis
gmx chi -s $PATHdata/$gro.gro -f $PATHdata/$trj.xtc -phi -psi -omega -all -maxchi 6
