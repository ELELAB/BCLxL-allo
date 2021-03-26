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
import os 

################################################
##########  FILTER CSVs          ###############
################################################
dataset = '18793'

os.makedirs('filtered_csv',exist_ok=True)

replicate1= '../../../../../replicate1/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_{}/csv_summary/summary_PPM.csv'.format(dataset)
replicate2= '../../../../../replicate2/CHARMM22star_TIP3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_{}/csv_summary/summary_PPM.csv'.format(dataset)
replicate31 = '../../../../../replicate3/CHARMM22star_TIPS3P/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_{}/csv_summary/summary_PPM.csv'.format(dataset)
replicate32 = '../../../../../replicate3/a99SB-disp/md/cs_pred/ppm/ppm_exp_analysis/1000_ns_{}/csv_summary/summary_PPM.csv'.format(dataset)

# list_replicates 
replicates = {replicate1:'replicate1', replicate2:'replicate2', replicate31:'replicate31', replicate32:'replicate32'}

for path, folder_name in replicates.items(): 
    # load data
    data = pd.read_csv(path, index_col = 0)

    # list of atom types 
    atomtypes = ['CA','C','CB','H','HA','N']
    
    # filtered atomtypes 
    data_backbone = data[data['atomtype'].isin(atomtypes)][['idPDB','atomtype','ChiSq']]
    
    # List of residues to filter 
    residues_high_CA = [55,104,59,129,119]
    residues_high_CB= [103,73,65,75,118]
    
    
    # data_filtered
    data_backbone.loc[data_backbone['idPDB'].isin(residues_high_CA), 'ChiSq'] = np.nan
    data_backbone.loc[data_backbone['idPDB'].isin(residues_high_CB), 'ChiSq'] = np.nan
    
     
    # aggregated data 
    summary = data_backbone.groupby('atomtype').agg('mean')['ChiSq']
    print(summary)
   # summary.reset_index(inplace=True)
    #summary['atomtype']=summary.index
    os.makedirs(f'./filtered_csv/{folder_name}', exist_ok=True)
    summary.to_csv(f'./filtered_csv/{folder_name}/out_chi_squared_BB_filtered.csv', index=True)


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
    replicate1 = './filtered_csv/replicate1/out_chi_squared_BB_filtered.csv'
    replicate2 =  './filtered_csv/replicate2/out_chi_squared_BB_filtered.csv'
    replicate31 =  './filtered_csv/replicate31/out_chi_squared_BB_filtered.csv'
    replicate32  = './filtered_csv/replicate32/out_chi_squared_BB_filtered.csv'

     # list_replicates 
    replicates = [replicate1, replicate2, replicate31, replicate32]
    
    # list replicates 
    data_replicates = [pd.read_csv(replicate).sort_values('atomtype') for replicate in replicates]
    data_replicates =[ df[df['atomtype'].isin(['CA','CB','C','N','H'])] for df in data_replicates]
    
    #labels 
    labels  = ['replicate 1 (CHARMM22* + TIP3P)',
            'replicate 2 (CHARMM22* + TIP3P)',
            'replicate 3 (CHARMM22* + TIPS3P)',
            'replicate 3 (a99SB-disp)']

       
    #labels bar
    atomtypes = data_replicates[0]['atomtype']
    labels_dic ={'N':'N', 'H':'H','CA':r'C$\mathregular{\alpha}$','CB':r'C$\mathregular{\beta}$', 'HA':r'H$\mathregular{\alpha}$','C':'C' }
    labels_bar = [labels_dic[atomtype] for atomtype in atomtypes]
    
    # load data. Creates a table where every row is an atom type (CA, CB ..) and every column is the agreggated red chi-sq for one trajectory. 
    # columns are named as in lables
    compilation_data = pd.DataFrame(data={labels[i]: data_replicates[i]['ChiSq'] for i in range(len(replicates))})
   
     
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
    ax.legend(prop=font_properties_small_text)
    fig.savefig('comparison_replicates_BB_{}_filtered.pdf'.format(dataset))
