This folder contains the files and the scripts to create heatmap plots with the local (per residue) reduced 
chi-squared (redChiSq) for all backbome atom types and for comparison with both chemical shift data set 2M04. 

The plots are made using a compilation with the backbone reduced redChiSq. The compilation was done manually using the output from 
ppm_exp_analysis stored in the following paths: 

>>../../../../../replicate1/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_18793/csv_summary/summary_PPM.csv
>>../../../../../replicate2/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_18793/csv_summary/summary_PPM.csv
>>../../../../../replicate3/CHARMM22star_TIPS3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_18793/csv_summary/summary_PPM.csv
>>../../../../../replicate3/a99SB-disp/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_18793/csv_summary/summary_PPM.csv



Plots are collected under folder plots. 

###########################################################################################
####                          REPRODUCE                                 ##################
##########################################################################################
0. Source Python 3.7 environment.
>>>  source /usr/local/envs/py37/bin/activate     

1. Execute:
>>> python plots_local_BB.py
