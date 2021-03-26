

import os
import biopandas
from biopandas.pdb import *
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

filename = 'torsions_all'
dataset = '18792'

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

# xvg file parser. From the output file of gmx sasa (resarea.xvg) creates a pandas dataframe of the form [residue | sasa | std]

import operator
import re
import numpy
import pandas

# A regular expression to split a string around spaces, while keeping quoted
# strings together. The regexp is modified from
# <http://stackoverflow.com/questions/2785755/how-to-split-but-ignore-separators-in-quoted-strings-in-python>
RE_SPLIT = re.compile(r'''((?:[^\s"']|"[^"]*"|'[^']*')+)''')


class XvgSingle(object):
    """
    A dataset with metadata from a XVG file
    A XVG file contains one or multiple dataset along with metadata and plot
    instructions. This class gives access to one dataset and its metadata.
    The accessible metadata are the dataset title, the label of the X axis, the
    label of the Y, and the name of the columns. They are accessible *via* the
    `title`, `xlabel`, `ylabel`, and `columns attributes, respectively.
    The data can be accessed and manipulated like a numpy array. In addition,
    columns can be acessed by their name.
    """
    _translate_commands = {
        'title': 'title',
        'xaxis': 'xlabel',
        'yaxis': 'ylabel'
    }

    def __init__(self):
        self._columns = {}
        self._array = np.empty((0,))
        self.xlabel = ''
        self.ylabel = ''
        self.title = ''

    def parse(self, content):
        """
        Read an input to fill the current instance
        Only the first dataset is read from the input. The reading stops at
        the first occurence of '//' in the file or at the end of the file if
        there is only one dataset.
        Parameters
        ----------
        content: iterator
            Iterator over the lines of a XVG file.
        """
        values = []
        for line in content:
            if line.startswith('//'):
                break  # For now we read only one dataset
            elif line.startswith('@'):
                self._parse_header_line(line)
            elif line.startswith('#'):
                continue
            else:
                values.append(line.split())
        self._array = numpy.array(values, dtype='float')

    @property
    def columns(self):
        items = list(sorted(self._columns.items(), key=operator.itemgetter(1)))
        columns = []
        i = 0
        while items:
            if items[0][1] == i:
                columns.append(items.pop(0)[0])
            else:
                columns.append('')
            i += 1
        return columns

    @classmethod
    def from_iter(cls, content):
        """
        Create an instance from an iterator over XVG lines
        """
        xvg = cls()
        xvg.parse(content)
        return xvg

    @classmethod
    def from_fname(cls, fname):
        """
        Create an instance from the path to a XVG file
        """
        xvg = cls()
        with open(fname, 'rt') as infile:
            xvg.parse(infile)
        return xvg

    def _parse_header_line(self, line):
        tokens = RE_SPLIT.findall(line[1:].strip())
        # the quotation symbols are kept, but we want them removed
        tokens = [
            token[1:-1]
            if token[0] == token[-1] and token[0] in '"\'' else token
            for token in tokens
        ]
        command = tokens.pop(0)
        if command == 'title':
            self.title = tokens.pop(0)
        elif command in self._translate_commands:
            setattr(self, self._translate_commands[command], tokens.pop(0))
        elif command.startswith('s') and tokens[0] == 'legend':
            tokens.pop(0)
            self._columns[tokens.pop()] = int(command[1:]) + 1

    def __getitem__(self, key):
        try:
            hash(key)
        except TypeError:
            key_is_mutable = False
        else:
            key_is_mutable = True
        if key_is_mutable and key in self._columns:
            return self._array[:, self._columns[key]]
        return self._array[key]

    def __getattr__(self, attr):
        return getattr(self._array, attr)


#content = open('resarea.xvg', 'r')
#content = content.readlines()
#
#
#file = XvgSingle()
#file.parse(content)
#
#array = file._array
#array[:,0] = array[:,0].astype(int)
#SASA = pd.DataFrame(data=array, columns=['residue','sasa','std'])


# creates auxiliary file using the function file_processor
file_processor(filename)

# creates pd.DataFrame with the population of vals and leus
pop = pd.read_csv('populations_{}.csv'.format(dataset), sep=';', header=None)[[0,1,2,3,4,5,6]]
pop = pop.rename(columns={0:'res_nr',
    1:'res_name',
    2: 'C1',
    3:'C2',
    4:'g_plus',
    5:'trans',
    6:'g_minus'})

print(pop)
pop['g_plus_sim'] = 0
pop['trans_sim']= 0 
pop['g_minus_sim'] = 0

print(pop)
# creates a pandas dataframe with the structure: time| torsion1 , torsion2, ... torsion n
data = pd.read_csv('file_modified',sep=' ', header=0, skipinitialspace=True)
#data = data.drop(columns=['t_637']) #I drop a column in the bh3
# extracts the column names except time
columns = data.columns[1:]
print(columns)
print(len(columns))
# extracst list of residues in the structure
path_ref= os.path.join(os.getcwd(),'reference.pdb')
newPDB=PandasPdb().read_pdb(path_ref) #read reference file
list_residues = newPDB.amino3to1()['residue_name'].to_list()
print(list_residues)
#atoms = newPDB.df['ATOM'] #pandas data frame with structure

# creates a dictionary with colnames and label
angles = {'V':r'$\chi ^1$','L':r'$\chi ^2$'}

# list of residues 
residues = [list_residues[int(columns[i].split('_')[-1])] for i in range(len(columns))]
residue_numbers = [int(columns[i].split('_')[-1]) for i in range(len(columns))]

new_columns = {str(columns[i]):'{} {}'.format(residue_numbers[i],residues[i])  for i in range(len(columns))}
data = data.rename(columns=new_columns)
print(new_columns)
print(data.columns)

#creates a dictionary with the SASA of the residues 
#SASA = {101:r'$ 0.045 \pm 0.027$',
#        86: r'$ 0.033 \pm 0.00 $',
#        115:  r'$ 0.590 \pm 0.372 $',
#        10:  r'$ 0.409 \pm 0.225 $',
#        95:  r'$ 0.170 \pm 0.241 $',
#        9:  r'$ 0 \pm 0 $',
#        72:  r'$ 0.164 \pm 0.081 $',
#        8: r'$ 0.085 \pm 0.030 $',
#        90: r'$ 0.013 \pm 0.019 $',
#        110: r'$ 0.210 \pm 0.110 $',
#        138: r'$ 0.027 \pm 0.038 $',
#        50: r'$ 0.049 \pm 0.032 $',
#        122: r'$ 0.053 \pm 0.038 $'}
#
# plots histograms
#residue_numbers.remove(134)
#nrplot = []
#for i in range(len(residue_numbers)):
#
#    print(data[data.columns[i+1]])
#    print(residue_numbers[i])
#    nr = i//8
#    nrplot.append(nr)
#    if i!=0:
#        if nrplot[-1] != nrplot[-2]:
#            plt.savefig('test_{}.pdf'.format(nr))
#    if i%8==0:
#        fig, axs = plt.subplots(4,2,figsize=(20,18))
#        axes = axs.flatten()
#           # axes[i%8].title.set_text('teeeeest') 
#    #axes[i%8].set_title
#
#    #plots histograms
#    data.hist(column = data.columns[i+1],
#            bins=36,
#            ax=axes[i%8],
#            alpha=0.5)
#    
##    data[data.columns[i+1]].plot.kde(ax=axes[i%8],secondary_y=True,sharex=False,legend=False)
#    axes[i%8].set_ylabel('frequency' + '||' + r'SASA $ {}\pm{} nm^2$'.format(float(SASA[SASA.residue == residue_numbers[i]]['sasa']),
#                float(SASA[SASA.residue == residue_numbers[i]]['std'])))
#    axes[i%8].set_xlabel(angles[residues[i]] + ' [rad]')
#    axes[i%8].set_yticks(np.arange(0,len(data[data.columns[i+1]])/2,10000))
#    axes[i%8].set_yticklabels([str(round(n,2)) for n in np.arange(0,len(data[data.columns[i+1]])/2,10000)/len(data[data.columns[i+1]])])
#    #add text
#    ax2 = axes[i%8].twinx()  # instantiate a second axes that shares the same x-axis instantiate a second axes that shares the same x-axis
#    ax2.plot(np.pi*np.array([63,174,-63])/180,[float(pop[pop[0] ==residue_numbers[i]][j]) for j in [4,5,6]],'o')
#    ax2.set_ylabel('Population prediction (from NMR)')
#    plt.savefig('test_{}.pdf'.format(nr+1))

for column in new_columns.values(): 
    res_nr = int(column.split(' ')[0])
    state_trans = data[(data[column].between(-np.pi,-2*np.pi/3)) | (data[column].between(2*np.pi/3,np.pi))][column]
    pop_trans = len(state_trans) / len(data[column])
    pop.at[pop[pop.res_nr == res_nr].index,'trans_sim'] = pop_trans

    state_g_plus = data[data[column].between(0,np.pi*2/3)][column]
    pop_g_plus = len(state_g_plus) / len(data[column])
    pop.at[pop[pop.res_nr == res_nr].index,'g_plus_sim'] = pop_g_plus

    state_g_minus = data[data[column].between(-np.pi*2/3,0)][column]
    pop_g_minus = len(state_g_minus) / len(data[column])
    pop.at[pop[pop.res_nr == res_nr].index,'g_minus_sim'] = pop_g_minus


    print(pop)
    print(column)
    print(pop_trans)
bad_res=[72,8,90,110,138,50,122,101,86,115,10,95,9]
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
        plt.savefig(res + '_' + 'a99sb_' +  dataset +  '.pdf')

pop.to_csv('populations_table_{}'.format(dataset))
    


