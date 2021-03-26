import pandas as pd 
import os
import numpy as np

#experimental 
path_exp='./experimental.csv'


#creates a pd.DataFrama to collect the input to Sider
data_sider = pd.DataFrame(columns=['residue_number','residue','C1','C2'])
print(data_sider)

data = pd.read_csv(path_exp)

#filter data to get only ILE, LEU and VAL entries.
data = data[data.residue.isin(['ILE','LEU','VAL'])][['idPDB','residue','atomtype','CS_exp']]

#drop residues with nan in idPDB
data = data.dropna(subset=['idPDB'])

#iterates over all residues in filtered data
for i in set(data['idPDB'].to_list()):

    #retrieve residue
    residue = data[data.idPDB==i]
   
    #list of atoms for which chemical shift entries are available
    atoms = set(residue['atomtype'].to_list())

    #Add to the dataframe data_sider only VAL and LEU entries for with both methyl C cs are available.
    if ((set(residue['residue'].to_list()) == {'VAL'}) & ({'CG1','CG2'}.issubset(atoms))):
        CG1_cs = float(residue[residue.atomtype == 'CG1']['CS_exp'])
        CG2_cs = float(residue[residue.atomtype == 'CG2']['CS_exp'])
        new_entry=pd.DataFrame(data =[[int(i),'VAL',CG1_cs,CG2_cs]], columns=['residue_number','residue','C1','C2'])
        data_sider  = data_sider.append(new_entry, ignore_index=True)
        
    if ((set(residue['residue'].to_list()) == {'LEU'}) & ({'CD1','CD2'}.issubset(atoms))):
        CD1_cs = float(residue[residue.atomtype == 'CD1']['CS_exp'])
        CD2_cs = float(residue[residue.atomtype == 'CD2']['CS_exp'])
        new_entry=pd.DataFrame(data =[[int(i),'LEU',CD1_cs,CD2_cs]], columns=['residue_number','residue','C1','C2'])
        data_sider  = data_sider.append(new_entry, ignore_index=True)
    
    #Add to the dataframe data_sider only ILE entries for with cs of CD is available
    if ((set(residue['residue'].to_list()) == {'ILE'}) & ({'CD1'}.issubset(atoms))):
        CD_cs = float(residue[residue.atomtype == 'CD1']['CS_exp'])
        new_entry=pd.DataFrame(data =[[int(i),'ILE',CD_cs,np.nan]], columns=['residue_number','residue','C1','C2'])
        data_sider  = data_sider.append(new_entry, ignore_index=True)
          
#Dump a txt file with the input for sideR
data_sider.to_csv('input_file_sider.txt', sep=' ', header=False, index=False)
#print(data_sider)
