# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 17:40:29 2020

@author: dioni
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 16:52:47 2020

@author: dioni
"""

import pandas as pd 
import biopandas
import os
import biopandas
from biopandas.pdb import *
import itertools 

#copy and rename template
os.popen('cp plumed_template.dat plumed_torsions.dat')

#open plumed.dat template
file_plumed = open('plumed_torsions.dat','a')

#print arg
print_arg = []

# extracting atoms for potential 
path_ref= os.path.join(os.getcwd(),'reference.pdb')
newPDB=PandasPdb().read_pdb(path_ref) #read reference file
atoms = newPDB.df['ATOM'] #pandas data frame with structure
leus_nr = list(set(atoms[atoms.residue_name == 'LEU']['residue_number'].to_list()))

# Residues to monitor
leus_nr = list(set(atoms[atoms.residue_name == 'LEU']['residue_number'].to_list()))
leus_atoms = ['CA','CB','CG','CD1']
file_plumed.write('\n####    LEU RESIDUES ####\n')
print(leus_nr)
for nr in leus_nr:
    leu  =  atoms[(atoms.residue_name == 'LEU') & (atoms.residue_number == nr)]
    print(leu)
    atoms_chi = [int(leu[leu.atom_name==j]['atom_number']) for j in leus_atoms]
    file_plumed.write('t_{}: TORSION ATOMS={},{},{},{}\n'.format(nr,atoms_chi[0],atoms_chi[1],atoms_chi[2],atoms_chi[3]))
    print_arg.append('t_{}'.format(nr))

# Residues to monitor
vals_nr = list(set(atoms[atoms.residue_name == 'VAL']['residue_number'].to_list()))
vals_atoms = ['N','CA','CB','CG1']
file_plumed.write('\n####    VAL RESIDUES ####\n')

for nr in vals_nr:
    val  =  atoms[(atoms.residue_name == 'VAL') & (atoms.residue_number == nr)]
    print(val)
    atoms_chi = [int(val[val.atom_name==j]['atom_number']) for j in vals_atoms]
    file_plumed.write('t_{}: TORSION ATOMS={},{},{},{}\n'.format(nr,atoms_chi[0],atoms_chi[1],atoms_chi[2],atoms_chi[3]))
    print_arg.append('t_{}'.format(nr))





# Residues to monitor
  
iles_nr = list(set(atoms[atoms.residue_name == 'ILE']['residue_number'].to_list()))
iles_atoms = ['CA','CB','CG1','CD']
file_plumed.write('\n#### ILE RESIDUES ####\n')

  
for nr in iles_nr:
    ile  =  atoms[(atoms.residue_name == 'ILE') & (atoms.residue_number == nr)]
    print(ile)
    atoms_chi = [int(ile[ile.atom_name==j]['atom_number']) for j in iles_atoms]
    file_plumed.write('t_{}: TORSION ATOMS={},{},{},{}\n'.format(nr,atoms_chi[0],atoms_chi[1],atoms_chi[2],atoms_chi[3]))
    print_arg.append('t_{}'.format(nr))


#PRINT ARG STATEMENT 
file_plumed.write('\n\nPRINT ARG=' + ','.join(print_arg) + ' FILE=torsion_all')

    

    








