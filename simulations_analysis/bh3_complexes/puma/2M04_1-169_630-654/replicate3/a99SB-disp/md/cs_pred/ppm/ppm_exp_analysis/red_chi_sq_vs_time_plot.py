# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 10:13:24 2019

@author: dioni
"""

import numpy as np
import sys
import re
import pandas as pd
import os 
import math
import re
import matplotlib.pyplot as plt

def read_output_ppm(filename):
    file=pd.read_csv(filename)
    data=file.values

    labels = file['atomtype']
    #print(labels)
    return data[:,2], labels


folders=np.arange(0,1001,step=25) #times of the subtrajectories: 0, 25, 50, ..., 1000 ns
print(folders)

#labels=['C','CB','CA','N','H','HA']
structures = [18793]

for struct in structures: 

    #Check the number of atoms for which output has been create since there might be missing bb atoms. 
    data0, labels=read_output_ppm('./0_ns_{}/csv_summary/output_chi_squared_BB.csv'.format(struct))
    num_atoms=np.size(data0, axis=0)

    #create empty array of dimension num of atom classes x number of datapoints. 
    data_all=np.zeros([num_atoms,len(folders)])

    for i in range(len(folders)):
        name_folder='{}_ns_{}'.format(folders[i],struct)
        os.chdir(name_folder)
        new_data, labels=read_output_ppm('./csv_summary/output_chi_squared_BB.csv')
        data_all[:,i]=new_data
        os.chdir('..')
    
    #create plots
    print(data_all)
    plt.figure()
    for j in range(len(labels)):
        plt.plot(folders,data_all[j,:],label=labels[j])
        plt.xlabel('time [ns]')
        plt.ylabel(r'$\chi ^2_{red}$ of chemical shift')
    
    plt.legend()
    plt.title(r'$\chi^2_{red}$ of backbone CS  predictions (free, a99SB-disp)') #change title of the plot
    plt.savefig('convergence_CS_BB_{}'.format(struct))
    
