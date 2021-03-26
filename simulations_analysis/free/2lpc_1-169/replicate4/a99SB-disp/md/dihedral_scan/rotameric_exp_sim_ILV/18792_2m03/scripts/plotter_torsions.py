

import os
import biopandas
from biopandas.pdb import *
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

filename = 'torsion_all'
dataset = '18792'
output_sider='output_sider.txt'


def file_processor(filename):
    '''function that processe the plumed output and creates an auxiliary file to be parsed to the read_csv 
    function of pandas'''

    file = open(filename,'r+')
    lines = file.readlines()
    newfile = open('file_modified','w')
    for i in range(len(lines)): 
        line = lines[i]
        line = line.replace('#! FIELDS', '')
        if i==0:
            newfile.write(line)
        else:
            if re.match(r'\s*\d+\.*\d*.*',line):
                newfile.write(line)
              #  print(line)    
    file.close()
    newfile.close()



def sider_output_parser(filename):
    '''This is a function that takes as an input a .csv file with the output of sider 
    [https://www.ucl.ac.uk/hansen-lab/sider_dfh.html] (with a copy of the whole page)
    and returns a .csv file of the form: resnr;residue_name;CS1;CS2;pop1;pop2;pop3)'''

    file = open(filename, 'r+')
    lines =file.readlines()
    newfile = open('populations_sider.csv','w')

    for i in range(len(lines)):
        line = lines[i]
        line = line.lstrip()
        line = line.rstrip()
        if re.match(r'\d+\s+.+',line):
            new_line = re.split(r'\s+',line)
            new_line = [element if element!='-' else '' for element in new_line]
            newfile.write(';'.join(new_line) + '\n')

    newfile.close()


# creates a auxiliary file from plumed output
file_processor(filename)

#creates the file populations_sider.csv
sider_output_parser(output_sider)


# creates pd.DataFrame with the populations estimated from sideR contained in the file populations_sider.csv
pop = pd.read_csv('populations_sider.csv', sep=';', header=None, names=['res_nr','res_name','C1','C2','g_plus','trans','g_minus']) #[[0,1,2,3,4,5,6]]
pop['g_plus_sim'] = 0
pop['trans_sim']= 0 
pop['g_minus_sim'] = 0


# creates a pandas dataframe with the structure: time| torsion1 , torsion2, ... torsion n
data = pd.read_csv('file_modified',sep=' ', header=0, skipinitialspace=True)

#data = data.drop(columns=['t_637']) #I drop a column in the bh3


# extracts the column names except time
columns = data.columns[1:]


# extracst list of residues in the structure
path_ref= os.path.join(os.getcwd(),'reference.pdb')
newPDB=PandasPdb().read_pdb(path_ref) #read reference file
list_residues = newPDB.amino3to1()['residue_name'].to_list()



# creates a dictionary with colnames and label
angles = {'V':r'$\chi ^1$','L':r'$\chi ^2$'}

# list of residues 
residues = [list_residues[int(columns[i].split('_')[-1])] for i in range(len(columns))]
residue_numbers = [int(columns[i].split('_')[-1]) for i in range(len(columns))]

#rename the columns of the dataframe to RESRN RES, e.g 101 V
new_columns = {str(columns[i]):'{} {}'.format(residue_numbers[i],residues[i])  for i in range(len(columns))}
data = data.rename(columns=new_columns)


# computes the rotameric populations from the simulation (from dataframe 'data') and fills out 
# the dataframe pop (containing the rotameric populations from sideR
for column in new_columns.values(): 
    res_nr = int(column.split(' ')[0]) #extracts the residue numer

    # trans : -180, -120 and 120, 180
    state_trans = data[(data[column].between(-np.pi,-2*np.pi/3)) | (data[column].between(2*np.pi/3,np.pi))][column]
    pop_trans = len(state_trans) / len(data[column])
    pop.at[pop[pop.res_nr == res_nr].index,'trans_sim'] = pop_trans

    # g plus = 0, 120
    state_g_plus = data[data[column].between(0,np.pi*2/3)][column]
    pop_g_plus = len(state_g_plus) / len(data[column])
    pop.at[pop[pop.res_nr == res_nr].index,'g_plus_sim'] = pop_g_plus

    # g_minus = -120, 0
    state_g_minus = data[data[column].between(-np.pi*2/3,0)][column]
    pop_g_minus = len(state_g_minus) / len(data[column])
    pop.at[pop[pop.res_nr == res_nr].index,'g_minus_sim'] = pop_g_minus



#bad_res=[72,8,90,110,138,50,122,101,86,115,10,95,9]


os.makedirs('../plots', exist_ok=True)


for res in ['val','leu','ile']:
    plt.figure()
    markers = ['o','v','*']
    states = ['g_plus','trans','g_minus']
    labels = [r'g$^+$','trans',r'g$^-$']
    states_sim = ['g_plus_sim','trans_sim','g_minus_sim']
    for i,j in zip(states,states_sim):
        
        x=pop[pop.res_name == res][i]
        y=pop[pop.res_name == res][j]
        tags = list(pop[pop.res_name == res]['res_nr'])
        plt.scatter(x,y,marker = markers[states.index(i)],label=labels[states.index(i)])
       # for l in range(len(tags)):
       #     if tags[l] in bad_res:
       #         plt.scatter(list(x)[l],list(y)[l],marker = markers[states.index(i)],color='red')

        plt.xlabel('Population from CS data', fontsize='x-large')
        plt.ylabel('Population from MD simulation',fontsize='x-large')
        plt.legend()
        for m, txt in enumerate(tags):
            if np.abs(list(x)[m]-list(y)[m]) > 0.20:
                plt.annotate(int(txt), (list(x)[m], list(y)[m]),fontsize='large')
        plt.plot(np.arange(0,1.1,0.05)+0.2,np.arange(0,1.1,0.05),linestyle='dotted',color='black')
        plt.plot(np.arange(0,1.1,0.05),np.arange(0,1.1,0.05),linestyle='--',color='black',linewidth=0.2)
        plt.plot(np.arange(0,1.1,0.05)-0.2,np.arange(0,1.1,0.05),linestyle='dotted',color='black')
        plt.xlim(-0.05,1.05)
        plt.ylim(-0.05,1.05)
        plt.title(res.upper())
        plt.savefig('../plots/' + res + '_' + dataset +  '.pdf')




pop.to_csv('table_comparison_exp_sim_{}.csv'.format(dataset))



#remove auxiliary files
os.system('rm file_modified populations_sider.csv')


