This file contains the python script to generate the plots comparing the backbone chemical shifts of three different chemical shift 
data sets used (2LPC [18250] and 2M03 [18792] for BCLxL free, and 2M04 for BCL-PUMA). 

It can be used to identify potential differences in chemical shift records for the same system coming from two different experiments or
two identify location where the binding of PUMA might induce changes at the level of the backbone in the BCL-xL.

To use the plots, we use the compilations csv files from three runs of the chemical shift comparative analysis, retrieve the experimental entries
for the data sets. 

Generation of the symbolic links: 

>>> ln -s ../../../simulations_analysis/free/2lpc_1-169/replicate1/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_18250/csv_summary/summary_PPM.csv summary_PPM_18250.csv

>>> ln -s ../../../simulations_analysis/free/2lpc_1-169/replicate1/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_18792/csv_summary/summary_PPM.csv summary_PPM_18792.csv

>>> ln -s ../../../simulations_analysis/bh3_complexes/puma/2M04_1-169_630-654/replicate1/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_18793/csv_summary/summary_PPM.csv summary_PPM_18793.csv


