


# source /usr/local/envs/py37/bin/activate 

#creates reference
path_traj=../../../../../../../../../../../simulations/bh3_bcl/bclxl/free/2lpc_1-169/replicate4/a99SB-disp/md/Mol_An/center_traj.xtc
path_tpr=../../../../../../../../../../../simulations/bh3_bcl/bclxl/free/2lpc_1-169/replicate4/a99SB-disp/md/md.tpr

echo "19\n"|gmx trjconv -f ${path_traj} -s ${path_tpr} -n index -b 0 -e 0 -o reference.pdb



for bmrbentry in 18250 18792 ; 
do
for i in 1000; # 0 25 50 75 100 125 150 175 200 225 250 275 300 325 350 375 400 425 450 475 500 525 550 575 600 625 650 675 700 725 750 775 800 825 850 875 900 925 950 975 1000 ;
do
	mkdir ${i}_ns_${bmrbentry}
	cd ${i}_ns_${bmrbentry}
	cp ../delta_CS.py .
	cp -r ../aux_files .
        pwd 

        python delta_CS.py -exp $bmrbentry  -ref ../reference.pdb -pdb_mapping yes -histograms no -ch3 ../../ch3shift_results/out_ch3shift_${i}.txt
	

	rm ./delta_CS.py 
#	rm ./output_*
	rm ./reference_*	
	rm -r ./aux_files
        cd ..
done
done
#-ch3 ..:/ch3shift/out_ch3shift_${i}.


rm ./#*

