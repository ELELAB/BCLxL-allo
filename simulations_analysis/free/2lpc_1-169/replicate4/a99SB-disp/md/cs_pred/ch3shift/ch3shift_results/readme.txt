#CH3shift for methyl chemical shift predictions. 

#Dependencies
#- CH3Shift_lib (folder with parameters from library) 
#- CH3Shift.R (R pipeline)
#- command.cmd (commands file)

#Input
# - traj_1000.pdb (ln -s  ../../../pdb_movies/traj_1000.pdb)

#To execute: 
>> tsp -N 1 R CMD BATCH CH3Shift.R 



#ALTERNATIVELY: 

#With the script ch3shift_collector.sh there is the possibility to collect the analysis for subtrajectories 
#(it will be very slow). To do this simply change the times (ns) accordingly in the for loop at the beginning
#of the file.
#ch3shift_collector.sh creates the index file from the centered full trajectory, 
#indicated at the top of the script (path_tpr, path_traj), and makes a subtrajectory with 1-ns spaced conformers out of it. 
#It modifies the file "command_template.cmd", a template of the command file that is changed for each subtrajectory. 
#The subtrajectories are created locally and deleted at the end of the analysis. 
#
# Once adjusted, simply run:

>> ch3shift_collector.sh




Last updated: 23.04.2020
