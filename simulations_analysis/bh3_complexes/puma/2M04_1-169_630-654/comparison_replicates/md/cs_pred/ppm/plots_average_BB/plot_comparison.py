#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 10:32:35 2020

@author: dionisio
"""


# libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
import seaborn as sns
from matplotlib.colors import ListedColormap
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
from matplotlib import rc

################################################
##########   GLOBAL PLOT SETTING ###############
################################################

# creates object to pass front properties to text in the plot 
path_to_helvetica ='./Helvetica.ttf'
font_properties = FontProperties(fname=path_to_helvetica,
        style='normal',
        stretch='normal')
font_properties_small_text = FontProperties(fname=path_to_helvetica,
        style='normal',
        stretch='normal',
        size='small')


# sets the font in math text to stix. Most "helvetica-looking" one from default options. 
rcParams['mathtext.fontset'] = 'stix'


 # color map settings (color blind)
cmap = ListedColormap(sns.color_palette("colorblind",8).as_hex())
idx_colors = [0,1,6,7,5,9,3,2,4,8] #indices from the blind 10 seaborn colorpalette used to match the color convention across figures
colors = sns.color_palette(np.array(sns.color_palette("colorblind",10))[idx_colors]).as_hex()




# set width of bar
barWidth = 0.08



#################################################
#########    DATA REFERENCES    #################
################################################

# cs_dataset
cs_datasets = ['18793']

# pdbs 
pdbs = ['2M04']



#################################################
#########    PLOTTING          #################
################################################


for j in range(len(cs_datasets)):
    dataset = cs_datasets[j]
    fig, ax = plt.subplots()

    # paths
    replicate1 = '../../../../../replicate1/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_BB.csv'.format(dataset)
    replicate2 = '../../../../../replicate2/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_BB.csv'.format(dataset)
    replicate31 = '../../../../../replicate3/CHARMM22star_TIPS3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_BB.csv'.format(dataset)
    replicate32 = '../../../../../replicate3/a99SB-disp/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_BB.csv'.format(dataset)

     # list_replicates 
    replicates = [replicate1, replicate2, replicate31, replicate32]
    
    # list replicate dataframes. We discard HA entries, since only two were available  
    data_replicates = [pd.read_csv(replicate).sort_values('atomtype') for replicate in replicates]
    data_replicates =[ df[df['atomtype'].isin(['CA','CB','C','N','H'])] for df in data_replicates]
        
    #labels 
    labels  = ['replicate 1 (CHARMM22* + TIP3P)',
            'replicate 2 (CHARMM22* + TIP3P)',
            'replicate 3 (CHARMM22* + TIPS3P)',
            'replicate 3 (a99SB-disp)']

       
    #labels bar
    atomtypes = data_replicates[0]['atomtype']
    print(atomtypes)
    labels_dic ={'N':'N', 'H':'H','CA':r'C$\mathregular{\alpha}$','CB':r'C$\mathregular{\beta}$', 'HA':r'H$\mathregular{\alpha}$','C':'C' }
    labels_bar = [labels_dic[atomtype] for atomtype in atomtypes]
    
    # load data. Creates a table where every row is an atom type (CA, CB ..) and every column is the agreggated red chi-sq for one trajectory. 
    # columns are named as in lables
    compilation_data = pd.DataFrame(data={labels[i]: data_replicates[i]['mean_red_ChiSq'] for i in range(len(replicates))})
   
     
   # plt.figure() #starts new figure 

    # Set position of bar on X axis
    r1 = np.arange(len(labels_bar))

    for i in range(len(replicates)):
        col =  compilation_data.columns[i] 
        ax.bar([x + barWidth*i for x in r1], list(compilation_data[col]), width=barWidth, edgecolor='white', label=labels[i], color=colors[i])
     
    # Labels, ticks, ticklabels and spine stettings 
    max_y = 7
    ax.set_xlabel('Atom type', fontproperties = font_properties)
    ax.set_ylabel(r'Average reduced $\mathregular{\chi^2}$', fontproperties = font_properties)
    ax.set_title(r'Comparison of average reduced {} of {} - PUMA replicates [{}]'.format('$\mathregular{\chi^2}$','Bcl-x$\mathregular{_L}$', pdbs[j]), fontproperties=font_properties ) 
    ax.set_xticks([r + barWidth*(len(replicates)//2) for r in range(len(labels_bar))])
    ax.set_yticks(range(max_y))
    ax.set_xticklabels(labels_bar, fontproperties=font_properties)
    ax.set_yticklabels(range(max_y), fontproperties = font_properties)
    ax.set_ylim((0,max_y))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Create legend & Show graphic
    ax.legend(loc=1, prop=font_properties_small_text)
    fig.savefig('comparison_replicates_BB_{}.pdf'.format(dataset))
