source /usr/local/envs/py37/bin/activate
#creates reference
path_traj=../../../../../../../../../simulations/free/2lpc_1-169/replicate1/CHARMM22star_TIPS3P/md/Mol_An/center_traj.xtc
path_tpr=../../../../../../../../../simulations/free/2lpc_1-169/replicate1/CHARMM22star_TIPS3P/md/md.tpr


echo "19\n"|gmx trjconv -f ${path_traj} -s ${path_tpr} -n index -b 0 -e 0 -o reference.pdb




for bmrbentry in 18250 18792 ; 
do
for i in 0 25 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 475 500 525 550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000 ; do
	mkdir ${i}_ns_${bmrbentry}
	cd ${i}_ns_${bmrbentry}
	cp ../delta_CS.py .
	cp -r ../aux_files .
        pwd 

	if [  ${i} = 1000 ] 
	then	
        python delta_CS.py -exp $bmrbentry -bp ../../ppm_results/ppm_${i}/bmrb_pre.dat  -ref ../reference.pdb -pdb_mapping yes -histograms yes
	else
	echo "not 1000"
        python delta_CS.py -exp $bmrbentry -bp ../../ppm_results/ppm_${i}/bmrb_pre.dat  -ref ../reference.pdb -pdb_mapping no -histograms no
	fi

	rm ./CS_compare.py 
	rm ./output_*
	rm ./reference_*	
	rm -r ./aux_files
        cd ..
done
done
#-ch3 ..:/ch3shift/out_ch3shift_${i}.

rm ./#*
