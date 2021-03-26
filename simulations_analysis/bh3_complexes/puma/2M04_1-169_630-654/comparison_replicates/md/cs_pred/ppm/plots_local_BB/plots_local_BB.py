import pandas as pd 
import matplotlib.pyplot as plt 
import seaborn as sns
import numpy as np 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from matplotlib.colors import ListedColormap
from matplotlib import rcParams
from matplotlib.font_manager import FontProperties
from matplotlib import rc


#############################################
####      PLOTTING FORMATING ################
#############################################
# creates object to pass front properties to text in the plot 
path_to_helvetica ='./Helvetica.ttf'
font_properties = FontProperties(fname=path_to_helvetica,
        style='normal',
        stretch='normal',
        size = 'large')

font_properties_small_text = FontProperties(fname=path_to_helvetica,
        style='normal',
        stretch='normal',
        size='small')


# sets the font in math text to stix. Most "helvetica-looking" one from default options. 
rcParams['mathtext.fontset'] = 'stix'

# Data 
data2m04= pd.read_csv('bclxl_free_2m04_compilation_bb_red_chisq.csv')

# labels to colums 
columns = {'Replicate1': 'Rep 1 (CHARMM22* + TIP3P)',
           'Replicate2': 'Rep 2 (CHARMM22* + TIP3P)' ,
           'Replicate3_c22star': 'Rep 3 (CHARMM22* + TIPS3P)',
           'Replicate3_a99sb_disp':'Rep 3 (a99SB-disp)'}

atoms={'CA': r'C$\mathregular{\alpha}$',
       'CB': r'C$\mathregular{\beta}$', 
       'C': r'C',
       'N': 'N',
       'HA': r'H$\mathregular{\alpha}$', 
       'H': r'H'}


def plotter_chis_bb(data, nameplot):
    os.makedirs('plots', exist_ok=True)
    max_sat = 3
    for i in range(len(atoms.keys())):
        fig, ax = plt.subplots(ncols=1, figsize=(18,4))
        atom = list(atoms.keys())[i]
        title = atom
        data_plot = []

        #ax = axes[i]
        for column in columns.keys():
            data_atoms = data[data.atomtype==atom]
            idPDB = data_atoms['idPDB'].astype(int)
            chis = data_atoms[column].astype(float) 
            chis[chis>max_sat]=max_sat
            new = [chis.to_list().copy() for i in range(max_sat)]
            data_plot +=new
        
        #makes array from the data from the different atomtypes
        data_plot = np.array(data_plot)        
        hm = ax.imshow(data_plot, cmap='jet')
        ax.set_title(atoms[atom], fontproperties = font_properties)

        # ticks 
        ax.set_yticks([1+max_sat*i for i in range(len(columns.values()))])
        ax.set_yticklabels(labels = columns.values(),fontproperties = font_properties, rotation=20)
        xticks_locations = [1,10,20,30,40,44,45,50,60,70,80,90,100,110,120,130,140,150,160,169]
        xticks_labels = [1,10,20,30,40,44,'      85 ',90,100,110,120,130,140,150,160,170,180,190,200,209]
        ax.set_xticks(xticks_locations)
        ax.set_xticklabels(labels = xticks_labels, fontproperties = font_properties)
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.gcf().subplots_adjust(left=0.15)
        bar = fig.colorbar(hm, cax=cax, ticks = [1,2,3])
        bar.ax.set_yticklabels([' 1', ' 2', '>3'])
        plt.savefig(os.path.join('plots','{}_{}.pdf'.format(nameplot,atom)))
        plt.close()



plotter_chis_bb(data2m04,'2m04' )

