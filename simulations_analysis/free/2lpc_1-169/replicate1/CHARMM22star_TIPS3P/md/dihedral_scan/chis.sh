source /usr/local/gromacs-5.1.2_plumed-2.3b/bin/GMXRC.bash
PATHdata='../../../../../../../../simulations/free/2lpc_1-169/replicate4/a99SB-disp/md/Mol_An2'
gro=update
trj=center_traj
mkdir chis
cd chis
gmx_mpi chi -s $PATHdata/$gro.gro -f $PATHdata/$trj.xtc -phi -psi -omega -all -maxchi 6
