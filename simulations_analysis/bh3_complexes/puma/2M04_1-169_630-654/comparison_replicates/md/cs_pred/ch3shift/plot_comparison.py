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
        weight='light',
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

    fig, ax = plt.subplots()
    dataset = cs_datasets[j]

    # paths
    replicate1 = '../../../../replicate1/CHARMM22star_TIP3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'.format(dataset)
    replicate2 = '../../../../replicate2/CHARMM22star_TIP3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'.format(dataset)
    replicate31 = '../../../../replicate3/CHARMM22star_TIPS3P/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'.format(dataset)
    replicate32 = '../../../../replicate3/a99SB-disp/md/cs_pred/ch3shift/ch3shift_exp_analysis/1000_ns_{}/csv_summary/output_chi_squared_CH3.csv'.format(dataset)

     # list_replicates 
    replicates = [replicate1, replicate2, replicate31, replicate32]
    data_replicates = [pd.read_csv(replicate).sort_values(['residue','atomtype']) for replicate in replicates]
    print(data_replicates)

    #labels 
    labels  = ['replicate 1 (CHARMM22* + TIP3P)',
            'replicate 2 (CHARMM22* + TIP3P)',
            'replicate 3 (CHARMM22* + TIPS3P)',
            'replicate 3 (a99SB-disp)']
    
    #labels bar
    labels_bar = [(row['residue'], row['atomtype']) for idx,row in data_replicates[0].iterrows()]
    
    labels_tag = {'CD':r'C$\mathregular{\delta}$',
            'CG2':r'C$\mathregular{\gamma_2}$',
            'CD1':r'C$\mathregular{\delta_1}$',
            'CD2':r'C$\mathregular{\delta_2}$',
            'CG1':r'C$\mathregular{\gamma_1}$'}

    # load data
    print(len(replicates),len(labels))
    compilation_data = pd.DataFrame(data={labels[i]: data_replicates[i]['mean_red_ChiSq'] for i in range(len(replicates))})


    # Set position of bar on X axis
    r1 = np.arange(len(labels_bar))

    for i in range(len(replicates)):
        col =  compilation_data.columns[i] 
        ax.bar([x + barWidth*i for x in r1], list(compilation_data[col]), width=barWidth, edgecolor='white', label=labels[i], color=colors[i])
     
    # Labels, ticks, ticklabels and spine stettings 
    max_y = 5
    ax.set_xlabel('Atom type', fontproperties = font_properties)
    ax.set_ylabel(r'Average reduced $\mathregular{\chi^2}$', fontproperties = font_properties)
    ax.set_title(r'Comparison of average reduced {} of {} - PUMA replicates [{}]'.format('$\mathregular{\chi^2}$','Bcl-x$\mathregular{_L}$', pdbs[j]), fontproperties=font_properties ) 
    ax.set_xticks([r + barWidth*(len(replicates)//2) for r in range(len(labels_bar))])
    ax.set_yticks(range(max_y))
    ax.set_xticklabels(['({}, {})'.format(label[0],labels_tag[label[1]]) for label in labels_bar],
            fontproperties=font_properties,
            rotation=45)
    ax.set_yticklabels(range(max_y), fontproperties = font_properties)
    ax.set_ylim((0,max_y))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    # Create legend & Show graphic
    ax.legend(prop=font_properties_small_text)
    fig.savefig('comparison_replicates_CH3_{}.pdf'.format(dataset), bbox_inches='tight')
