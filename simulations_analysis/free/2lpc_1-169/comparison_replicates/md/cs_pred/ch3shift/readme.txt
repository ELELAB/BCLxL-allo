This folder contains the plots comparing the results of the average reduced chi squared for methyl C atoms computed from 
the trajectories of the different replicates with the goal of assessing the degree of agreement of each model with 
the available chemical shift data. The relative paths from where the data are fetched are included in the python script plot_comparison.py
and they are (placeholder for 18792 or 18250 depending on the data set used for comparison)


'../../../../replicate1/CHARMM22star_TIP3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'
'../../../../replicate1/CHARMM22star_TIPS3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv
'../../../../replicate2/CHARMM22star_TIP3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'
'../../../../replicate2/CHARMM22star_TIPS3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'
../../../../replicate3/CHARMM22star_TIP3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'
'../../../../replicate4/CHARMM22star_TIPS3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'
'../../../../replicate4/a99SB-disp/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'


###################################
####### REPRODUCE             #####
##################################
 0. source /usr/local/envs/py37/bin/activate
1. python plot_comparison.py 
