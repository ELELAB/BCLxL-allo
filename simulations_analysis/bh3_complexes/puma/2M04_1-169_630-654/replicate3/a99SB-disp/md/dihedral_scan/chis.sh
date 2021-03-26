source /usr/local/gromacs-5.1.2_plumed-2.3b/bin/GMXRC.bash


PATHdata='../../../../../../../../../simulations/bh3_complexes/puma/2M04_1-169_630-654/replicate1/CHARMM22star_TIP3P/md/Mol_An'
gro=update
trj=traj_centered
mkdir chis
cd chis
gmx chi -s $PATHdata/$gro.gro -f $PATHdata/$trj.xtc -phi -psi -omega -all -maxchi 6
