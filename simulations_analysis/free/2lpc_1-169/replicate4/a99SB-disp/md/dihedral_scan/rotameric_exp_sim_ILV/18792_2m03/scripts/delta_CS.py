# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 13:24:11 2019

@author: dioni
"""
"""
CS_delta: comparison of experimental chemical shift data with predictions from
PPM and ch3shift. 
    Copyright (C) 2019 Dionisio Sanchez Cabrera (dionisio.sanchez.cab@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
import numpy as np
import sys
import re
import gromacs
import pandas as pd
import os 
import math
import re
import Bio
from Bio import *
from Bio import pairwise2
from Bio import SeqUtils
import biopandas
from biopandas.pdb import *
import argparse
import pynmrstar
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

###############################################################################
###############              PARSING ARGUMENTS              ###################
###############################################################################

#### parsing argument
parser = argparse.ArgumentParser(description='Comparison of chemical shifts with experiments.')
parser.add_argument('-exp','--experimental',type=str, help='BMRB entry number of the entity OR name of the BMRB file.')
parser.add_argument('-bp','--bmrb_pre',type=str, help='Name of the bmrb-like file generated by ppmone.', default=None)
parser.add_argument('-ch3shift','--ch3_shift',type=str, help='Name of the output generated by ch3shift.',default=None)
parser.add_argument('-ref','--reference',type=str, help='Name of the pdb-file used as reference')

# --> not available. Currently, only charmm nomenclature accepted. Nomenclature of reference changed with pdb2gmx. 
# --> No capping goups accepted. Tested with a99SB-disp and CHARMM22* nomenclature in the provided reference.pdb
#parser.add_argument('-pdb_nom','--pdb_nom',type=int, help='Version of the PDB nomenclature'
#                    'used in your reference structure. Check the nomenclature.txt file and choose 1 if it is using the one under PDBV1 or 2 for the one under PDBV2. You can modify single entries to correct for discrepancies',default=2)
parser.add_argument('-pdb_mapping','--pdb_mapping',type=str, help='True if you want the reduced chi-squared results to be mapped onto a pdb. Default is yes', default='yes')
parser.add_argument('-histograms','--histograms',type=str, help='True if you want  histograms with backbone chi-squares. Default is yes', default='yes')
args=parser.parse_args()
parser.print_help()




###############################################################################
###############              FUNCTIONS DEFINITION              ################
###############################################################################

def process_bmrb(filename,name):
    try:
        file_pr=open(str(name),'w')
        file=open(filename,'r')
    except FileNotFoundError as err: 
        print('There is no ppm output file under that name')


    lines=file.readlines()
    start=0
    for line in lines: 
        if start==True:
            if re.match(r'\s+[0-9]+\s+[0-9]+.+',line):
                file_pr.write(line)
        elif re.match(r'\s+_Chem_shift_ambiguity_code\s*',line):
            start=True
    file_pr.close()
    
    
    
def extract_sequence(index):
    '''This function extracts the residue sequence from the a pandas dataframe where having 
    the columns 'resID' and 'residue'. '''
    col1=index.columns[0]
    col2=index.columns[1]
    
    resID_pre=index.iloc[0][col1]
    residue_pre=index.iloc[1][col2]
    
    list_residues_df=[residue_pre]
    
    for i, row in index.iterrows():
         if row[col1]!=resID_pre:
             list_residues_df.append(row[col2])
             resID_pre=row[col1]
    return Bio.SeqUtils.seq1(''.join(list_residues_df))


                
def renumber_sequentially(index,col1,col2):
    '''This function takes a pandas dataframe / index object and renumbers its (numerical)
    feature 'col' sequentially. E.g. for a sequence with gaps in the numbering it will create a new column 'id_seq' 
    with a sequential numbering starting at 0. This avoids troubles when to chains
    are present'''
    #obtain first residue and first residues id 
    resID_pre=index.iloc[0][col1]
    residue_pre=index.iloc[1][col2]
    
    
    list_ids=[] #list of the new values of col1 (i.e.=new residues_id)
    
    #Counter used to change residues upon iteration
    counter=0
    
    #List of the residues. Used to obtain the sequence. Second object returned by 
    #the function
    list_residues=[residue_pre]
    idx=0 #the new resID will start at 0
    for i, row in index.iterrows():
         if row[col1]!=resID_pre:
             list_ids=list_ids + counter*[idx]
             counter=1
             idx+=1
             list_residues.append(row[col2])
             resID_pre=row[col1]

         else:
             counter+=1
    list_ids=list_ids + counter*[idx]
    
    index.at[:,'id_seq']=list_ids

    return index



def experimental_df(entry):
    '''Takes a bmrb entry number or a NMR-star file and creates a pandas DataFrame'''

    if os.path.isfile(entry):
        print("You are submitting a .txt data file with experimental chemical shift entries")
        try:
            entryExp=pynmrstar.Entry.from_file(entry)
        except FileNotFoundError as err: #Hd
            print('The entry file you submitted is not valid.\n',err)
            sys.exit('The entry file you submitted is not a valid one. Please, check it that you are submitting the right file.')

    elif re.match('[0-9]+',entry):
        entry=int(entry)
        print("You are trying to retrieve the chemical shift dataset deposited in the BMRB data base under the number {}.".format(entry))
        try:
            entryExp= pynmrstar.Entry.from_database(entry)
        except OSError as err:
            print('The entry number you submitted is not right. ',err)
            sys.exit('There is no data set under the number alias you have provided. Check for mistakes and try again.')

    else: 
        raise Exception('You must provide a file path or a bmrb entry code')
        sys.exit('Provide a right input in the form of bmbr entry number or NMR-star file')

    #defines a list with all the atoms in the structure
    cs_results_sets=[]
    cs_results_sets_list=[]
    
    #iterates over all the atoms in Entryexp
    for chemical_shift_loop in entryExp.get_loops_by_category("Atom_chem_shift"):
        cs_results_sets.append(chemical_shift_loop.get_tag(['Comp_index_ID','Comp_ID','Atom_ID','Atom_type','Val','Val_err']))
        
    #for nmrstar entries of polymers. EntryExp will be a list with the 
    # structure [[part1],[part2]...] where each part is a list with all the atoms. 
    # This transforms the "list of lists" structure into a single list. 
    for m in range(len(cs_results_sets)):  
        for i in range(len(cs_results_sets[m])):
            line=cs_results_sets[m][i]
            cs_results_sets[m][i]=[float(j) if re.match(r'\d+\.*\d*',j) else j for j in line]
            cs_results_sets_list.append(line)
    
    #creates a dataframe and corrects datatypes. 
    df=pd.DataFrame(cs_results_sets_list,columns=['resID_old','residue','atomtype','element','CS_exp','error_exp']) 
    df['resID_old']=df.resID_old.astype(int)
    df['CS_exp']=df.CS_exp.astype(float)
    df['error_exp']=df.error_exp.astype(float)
    
    #creates the column id_seq, i.e. a sequential numbering of the residues
    df=renumber_sequentially(df,'resID_old','residue')
    sequence=extract_sequence(df.loc[:,['resID_old','residue']]) 
    return df,sequence
    


#def get_starts(index): 
#    '''This function takes an index object with the structure |idx|resID|residue| and returns an index object whose entries are the start of each residue. 
#    The names of the columns of the input index object does not matter. Important is only the order of the columns'''
#    #get columns names: 
#    col1=index.columns[0]
#    col2=index.columns[1]
#    
#    #retrieve starting resID and residue number 
#    resID_pre=index.iloc[0][col1]
#    residue_pre=index.iloc[1][col2]
#    
#    #creates empty index object/dataframe
#    df=pd.DataFrame([[resID_pre,residue_pre]],columns=[col1,col2])
#    
#    
#    for i, row in index.iterrows():
#         if row[col1]!=resID_pre:
#             print(row)
#             print('x')
#             df=df.append(row,ignore_index=False)
#             resID_pre=row[col1]
#    return df



def predicted_ppm_df(filename):
    '''predicted_df returns a dataframe with the predictions from ppm contained in the file 'filename' AND
    a the 1-letter sequence of the sequence or residues in the form: dataframe, seq'''
    
    path=os.path.abspath(filename)
    
    #processes output bmrb file from pp by trimming unncessary lines
    process_bmrb(path,'bmrb_pre_pro.dat') 

    #creates pandas dataframe and trims unnessary columns
    data=pd.read_csv('bmrb_pre_pro.dat', delim_whitespace=True,header=None)
    
    
    df=data.loc[:, [1,2,4,5]]    
    df.columns = ['resID_old','residue','atomtype','CS']
    df=renumber_sequentially(df,'resID_old','residue')
    sequence=extract_sequence(df.loc[:,['resID_old','residue']]) 

    return df,sequence



def predicted_ch3(filename):
    '''predicted_ch3 returns a dataframe with the predictions of ch3shift.'''
    #open the file
    path=os.path.abspath(filename)
    try:
        file=open(path, 'r')
    except FileNotFoundError as err: 
        print('The file with ch3shift that you have given does not exist') 

    #read lines
    lines=file.readlines()
    trimmed_lines=[] #list with the new modified lines
    
    
    #trimming the lines
    for line in lines:
        if re.search('\ANOTE',line)==None:
            newLine=line.replace('RESULT: ;','')
            newLine=newLine.rstrip()
            newLine=re.split(r'\s+;\s+',newLine)
            for i in range(len(newLine)):
                element=newLine[i]
                newLine[i]=element.lstrip()    
            trimmed_lines.append(newLine)
    df=pd.DataFrame(trimmed_lines) #pd dataframe from new lines
    df=df.loc[:, [0,1,2,4,5,6,7]] #trimmig unnecesary information

    df.columns=['residue','idPDB','chain','atomtype','CS','err1','err2'] #labeling 
    
    #changing data types and some changes when necessary.
    df.idPDB=df.idPDB.astype(int)
    df.CS[df.CS=='xxxxx']=np.nan
    df.CS=df.CS.astype(float)
    df.err1=df.err1.astype(float)
    df.err2=df.err2.astype(float)
    df.loc[df.chain=='UN', 'chain']= 'A'
    
    return df
    


def no_match_index(alignment):
    '''This function takes an alignment object and return the indices where there is
    discrepancy'''
    no_matches = []
    for i, (a, b) in enumerate(zip(alignment[0], alignment[1])):
        if a != b:
            no_matches.append(i)
    return no_matches



def match_index(alignment):
    '''This function takes an alignment object and return the indices where there is NO
    discrepancy'''
    matches = []
    for i, (a, b) in enumerate(zip(alignment[0], alignment[1])):
        if a == b:
            matches.append(i)
    return matches



def alignment(seq1,seq2,return_matches=False):
    '''alignment returns the consensus string of two sequences. On request (return_matches=True), it can also return
    the indices in the original list of residues that the residues in the consensus correspond to.'''
    
 
    alignments = pairwise2.align.globalxx(seq1,seq2)
    
    #Determine the positions of the differences
    m=no_match_index(alignments[0])
    
    #Retrieve the common sequence 'consensus'
    consensus=''
    resume=0
    for i in range(len(m)):
        consensus = consensus + alignments[0][0][resume:m[i]]
        resume=m[i]+1
    consensus=consensus + alignments[0][0][resume:]
    
    #What indices in my primitive predicted sequences does the consensus residues correspond to?
    alignment_seq2_original_align=pairwise2.align.globalxx(seq2,consensus)
    matches_seq2=match_index(alignment_seq2_original_align[0])
    
    #What indices in my primitive experimental sequences does the consensus residues correspond to?
    alignment_seq1_original_align=pairwise2.align.globalxx(seq1,consensus)
    matches_seq1= match_index(alignment_seq1_original_align[0])
    
   # assert(type(return_matches)==bool,'return_matches must be True or False')
    if return_matches==True:
        return consensus,matches_seq1,matches_seq2
    else:  
        return consensus
    
    
    
def check_nomenclature(reference,versionPDB='PDBV2'):
    '''check_nomenclature asserts that the nomenclature for the atom types used in your reference pdb matches
    the one in the nomenclature.txt file utilized thoughout the program for atom type conversion. 
    In case of mismatches the name of the atom and the residue will be promted for the user to change the 
    corresponding entry manually in the nomenculate.txt file''' 'ADD RETURN MISMATCHES'''
    
    table=pd.read_csv(os.path.join(os.getcwd(),'aux_files','nomenclature.txt'),delim_whitespace=True) #load conversion table
    residues_checked=[] #in oder not to check residues twice
    for idx,row in reference.iterrows():
        if row['residue_name'] not in residues_checked: #check if that residue has already been inspected. 
            atoms_list=reference[(reference.residue_name==row['residue_name']) & (reference.residue_number==row['residue_number'])]['atom_name'].tolist()
            for atom in atoms_list: 
                if atom not in table[table.residue==row['residue_name']][versionPDB].to_list():
                    print('The atom {} from residue {} is not part of the nomenclature you have given. Please, change the nomenclature file accordingly or choose another PDBversion.'.format(atom,row['residue_name']))
            residues_checked.append(row['residue_name'])
        
        

def bmrb_to_pdb(residue,atomtype,versionPDB='PDBV2'):
    
    '''bmbr_to_pdb takes the atomtype of residue in NMR-star nomenclature and translates
    it into pdb nomenculature. The user should check the version of the PDB''' 
    
    table=pd.read_csv(os.path.join(os.getcwd(),'aux_files','nomenclature.txt'),delim_whitespace=True)
    try:
        atomInPDB=table[(table.residue==residue) & (table.BMRB==atomtype)][versionPDB]
    except:
        print('BMRB atomtype not found in the nomenclature directory')

    return atomInPDB.to_list()[0]



def create_reference_right_nomenclature(ref, **kwargs):
    """Create a new reference PDB (reference_renamed.pdb) with using CHARMM27 nomenclature"""

        #calls pdb2gmx to rename atoms according to CHARMM27 nomenclature. It is relevant mostly for protons. 
    Pdb2gmx_ref=gromacs.tools.Pdb2gmx(
            f=ref,
            o='reference_renamed.pdb',
            ff='charmm27', #selects force field charmm27
            ignh=True, #ignore the hydrogens of the original pdb and insert new ones with CHARMM27 nomenclature
            input=("7"), #select no water model
            ) 
    Pdb2gmx_ref.run()
    
    #remove *.itp and *.top files from pdb2gmx
    pdb2gmx_output_to_remove = [f for f in os.listdir() if ((f.endswith('.itp')) or (f.endswith('.top')))]
    [os.remove(i) for i in pdb2gmx_output_to_remove]



def bmrb_to_pdb2(row):
    '''This function takes a row and returns the PDB nomenclature for that atomtype. It handles
    especial cases with equivalent protons. E.g. ppm might generate just one prediction for 
    Ala HB1 and Ala HB2 under the name Ala HB. This function corrects that. The PDB nomeclature is
    returns as a list, e.g.: ['CA'], ['HB1','HB2]'''
    
    table=pd.read_csv(os.path.join(os.getcwd(),'aux_files','nomenclature.txt'),delim_whitespace=True)
    try:
        
        atomInPDB=table[(table.residue==row['residue']) & (table.BMRB==row['atomtype'])][versionPDB]
        return [atomInPDB.to_list()[0]]

    except:
        try: 
            h_type = re.match(r'H[B-Z]+[1-9]*',row['atomtype'])[0]
            if h_type :
                matches_in_nom = table[table.BMRB.str.match(h_type) & (table.residue == row['residue'])]
                list_matches = matches_in_nom[versionPDB].to_list()
                return list_matches
        except:
            print('BMRB atom type {} not found in the nomenclature directory'.format(row['atomtype']))
            return None

#    print(atomInPDB)
    


def pdb_to_bmrb(residue,atomtype,versionPDB=3):
    
    '''bmbr_to_pdb takes the atomtype of residue in NMR-star nomenclature and translates
    it into pdb nomenculature. The user should check the version of the PDB''' 
    
    table=pd.read_csv(os.path.join(os.getcwd(),'aux_files','nomenclature.txt'),delim_whitespace=True)
    if versionPDB==3:
        try:
            atomInBMRB=table[(table.residue==residue) & (table.PDBV2==atomtype)]['BMRB']
        except:
            print('PDB atomtype not found in the nomenclature directory')
    elif versionPDB==2:
        try:
            atomInBMRB=table[(table.residue==residue) & (table.PDBV1==atomtype)]['BMRB']
        except:
            print('PDB atomtype not found in the nomenclature directory')
    print('The atomtype {} from residue {}  was changed to {}.'.format(atomtype,residue,atomInBMRB))
    return atomInBMRB.to_list()[0]
        


def renumbering_to_reference(reference,df):
    '''This function takes a reference PDB and a pandas Dataframe (experimental o predicted_ppm)
    and will create a new column names 'ResID_ref' that will assign each residue in the dataframe 
    the resID of the same residue in the reference pdf. This makes atoms/residue referencing easier across the code since
    ppm output and nmrstar file do not necessarily have the same residue numbering as our reference pdb. If not all
    residues are present in the pdb, the value of ResID_ref for that residue will be nan.
    Returns the modifed dataframe'''
    
    try:
        referencepdb=PandasPdb().read_pdb(reference)
        reference_df=referencepdb.df['ATOM']
        reference_df=renumber_sequentially(reference_df,'residue_number','residue_name')

    except FileNotFoundError as err:
        print(err, '\n', 'The path of the reference PDB is not right. Please, check again. ')
        
    #Extracts the sequences of the reference and the dataframe
    seq_ref=''.join(referencepdb.amino3to1()['residue_name']) #sequency reference pdb 
    seq_df=extract_sequence(df)
#    print(seq_ref,seq_df)

    #get consensus string of df and reference pdb 
    consensus,matches_df,matches_ref=alignment(seq_df,seq_ref,return_matches=True)
#    print('The consensus sequence is ',consensus)
#    print(matches_df,matches_ref)

    #renumbers
    for i in range(len(consensus)):
        id_ref=matches_ref[i]
        id_df=matches_df[i]
        idPDB=int(list(reference_df[reference_df.id_seq==matches_ref[i]]['residue_number'])[0])
        df.at[df.id_seq==matches_df[i],'idPDB']=int(idPDB)
        
    return df
        


def add_error_ppm(index):
    '''This function adds to 'copy' (the compiled_df) the error of the predictor to calculate
    the chi-square later on'''
    copy=index.copy()
    bb_atoms=['C','CB','CA','N','H','HA']
    listError={'C':1.44, 'CA':0.9, 'N':2.31, 'H': 0.43, 'HA':0.24, 'CB':1}
    
    
    error_protons=pd.read_csv(os.path.join(os.getcwd(),'aux_files','error_protons_ppm.txt'),names=['residue','Atomtype','ErrorPPM','Shifts','CH3shift/ArShift','ExperimentalRMSD'])
    error_protons.Atomtype=error_protons.Atomtype.str.replace('*','').tolist()
    for idx,row in copy.iterrows():
        if row['atomtype'] in bb_atoms:
            copy.at[idx,'error_predictor']=listError[row['atomtype']]
        elif re.match(r'H[B-Z]+[1-9]*',row['atomtype'] )!=None:
                if row['atomtype'] not in error_protons[error_protons.residue==row['residue']]['Atomtype'].tolist():
                    error_proton=error_protons[(error_protons.residue==row['residue']) & (error_protons.Atomtype==row['atomtype'][:-1])]['ErrorPPM']
                    
                elif row['atomtype'] in error_protons[error_protons.residue==row['residue']]['Atomtype'].tolist():
                    error_proton=error_protons[(error_protons.residue==row['residue']) & (error_protons.Atomtype==row['atomtype'])]['ErrorPPM']
                
                try:
                    copy.at[idx,'error_predictor']=error_proton
        
                except:
        #                print('No experimental data for this residue')
                    pass

    return copy

            
def eliminate_duplicate_hydrogens(index):
    '''This function should correct for the presence of equivalent redundant protons. E.g. HD1 and HD11,HD12,HD13 '''
    try:
        print('Eliminating duplicated hydrogens resulting from different nomenclatures ...')
        copy=index.copy()
        deleted=[]
        for idPDB in list(set(copy.idPDB)):
            res=copy[copy.idPDB==idPDB]
            unique_elements, counts_elements = np.unique(res['CS_exp'], return_counts=True)
            for i in range(len(unique_elements)):
                if counts_elements[i]>1:
    #                print(counts_elements[i])
                    #check
                    element=unique_elements[i]
                    atomtype=res[res.CS_exp==element]['atomtype'].to_list()[0] #first of the redundants H
                    atomtype_class=''.join([i for i in atomtype if not i.isdigit()])
    #                print(atomtype_class)
                    origin=res[res.CS_exp==element]['_merge'].to_list()[0] 
    #                print(origin)#origine when merging 
                    for idx,row in res[res._merge!=origin].iterrows(): 
    #                    print('test')
                        if re.match(atomtype_class,row['atomtype']):
    #                        print('this actually happens')
    #                        print(row)
                            if not row['CS_exp']>0: #is nan?
                                sel=res[res.CS_exp==element]
                                copy.at[sel.index,'CS']=row['CS']
                                print('eliminates {},{},{}.'.format(row['residue'],row['atomtype'],row['idPDB']))
                                
                                deleted.append(idx)

    #                            print(' CS_exp of {} in residue {} from {} to {}'.format(row['atomtype'],row['residue'],'nan',unique_elements[i]))
    #                            copy.at[idx,'CS_exp']=unique_elements[i]
        print('Done eliminating duplicated protons..')
        copy=copy.drop(deleted)
        copy=copy.reset_index(drop=True)
    except Exception as err:
        print('Something went trom when eliminating duplicated hydrogens')
        print(err)
                                
    return copy



def add_experimental_data_to_ch3shift(ch3shift,experimental):
    '''This function complements the database of the predictions of CH3shift 
    with the experimental data'''
    print('Adding experimental data to CH3shift predictions ...')
    ch3shift_copy=ch3shift.copy()
    for idx,row in ch3shift_copy.iterrows():
        
        idPDBue_exp=experimental[experimental.idPDB==row['idPDB']]
        sel=idPDBue_exp.atomtype.str.match(row['atomtype'])
        try:
            CS_exp=idPDBue_exp[sel]['CS_exp'].to_list()[0]
    
#            if not CS_exp >0.0: 
#                print('There is no experimental information for the prediction of the atom {} of the idPDBue {}.'.format(row['atomtype'],row['idPDBue']))
            ch3shift_copy.at[idx,'CS_exp']=CS_exp
        except IndexError:
            print('There is no experimental information for the prediction of the atom {} of the residue {} with number {}.'.format(row['atomtype'],row['residue'],row['idPDB']))
    print('Done with the addition of experimental data to CH3shift predictions.')
    return ch3shift_copy
                            
  
   
def histogram_plotter(data, step=0.1, feature='Chi-Sq reduced', name='histogram_plot',title='Histogram'): 
    data = np.array(data)
    bins_max = min((max(data),5))
    label = '{:.2f} % of predictions'.format(len(data[data<bins_max])/len(data) * 100)
    bins = np.arange(0, bins_max, step)
    data_entries, bins = np.histogram(data, bins = bins)
    binscenters = np.array( [0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)] ) 
    plt.figure()
    plt.bar(binscenters,data_entries/sum(data_entries),width= 0.8 * ( (bins_max)/len(binscenters) ), color = 'navy',label=label)
    plt.xlabel(feature)
    plt.ylabel('frequency')
    plt.legend()
    plt.title(title)
    plt.savefig(name)
    plt.close()
    
    
    

def plotter(df, reference, atomtype):
    reference_copy = reference.copy()
    reference_copy['b_factor'] = -1
    df_selection = df[(df.atomtype == atomtype) & (df.ChiSq>0)]
#    print(df_selection.columns)

    for idx,row in  df_selection.iterrows():
        try:
            idPDB, resi= row['idPDB'], row['residue']
#            print(idPDB,row['atom_name'][0],residue)
            sel = reference_copy[(reference_copy.residue_number == idPDB) & (reference_copy.residue_name == resi)]
            reference_copy.at[sel.index,'b_factor'] = row['ChiSq']
#            print(idPDB,resi,row['ChiSq'],sel.index)
        except: 
            print('I could not find in the PDB the atom {} of residue {}{} when mapping the ChiSq'.format(row['atomtype'], row['idPDB'], row['residue']))
    return reference_copy
        #set the b_factor in reference to that value
                
        

class visualizer ():
    def __init__(self, reference=None, ppm_data=None , ch3_data=None): 
        self.ppm_data = ppm_data
        self.ch3_data = ch3_data
        self.reference = reference

    def PDBmapping(self, plotarg='BB'):
        '''maps results onto a pdb for visualization. Depends on function "plotter" '''
        
        if plotarg=='BB':
            df = self.ppm_data
            df= df[df.atomtype.isin(['C','CA','CB','HA','H','N']) & (df.ChiSq > 0.000)]
            os.makedirs('BB_pdbs', exist_ok=True)

            for atom in ['C','CA','CB','HA','H','N']:
                new_reference = plotter(df, self.reference, atom)
                newPDB = PandasPdb().read_pdb(args.reference)
                newPDB.df['ATOM']=new_reference 
                path = os.path.join('.','BB_pdbs', '{}_chi_squared_ppm.pdb'.format(atom))
                print(path)
                newPDB.to_pdb(path = path) 
        
                
#
        elif plotarg=='CH3_carbons':
            carbons_ch3 = {'LEU':['CD1','CD2'],'ILE':['CG2','CD'],'VAL':['CG1','CG2'],'THR':['CG2']}
            assert self.ch3_data is not None, 'There is no data from CH3shift'
            df = self.ch3_data
            os.makedirs('CH3_pdbs', exist_ok=True)
            for residue in carbons_ch3.keys():
                for C in carbons_ch3[residue]:
                    df_sel = df[(df.residue == residue) & (df.atomtype == C)] #only the carbon of interest 
                    new_reference = plotter(df_sel,self.reference,C)
                    newPDB = PandasPdb().read_pdb(args.reference)
                    newPDB.df['ATOM']=new_reference 
                    path = os.path.join('.','CH3_pdbs', '{}_{}_chi_squared_ppm.pdb'.format(residue,C))
                    newPDB.to_pdb(path = path)
                    
                    
    def histogramer_plotter(self):
        os.makedirs('histograms_BB', exist_ok=True)
        for atom in ['C','CA','CB','HA','H','N']:
            df = self.ppm_data
            df = df[(df.atomtype==atom) & (df.ChiSq>0.000)]
            try:
                histogram_plotter(df['ChiSq'], step=0.1, feature=r'$\chi^2_{red}$', name='./histograms_BB/histogram_plot_{}'.format(atom),title='{}'.format(atom))
            except: 
                print('There was a problem with the histogram for the atom type {}. Maybe there are no data for it.'.format(atom))

                        
                        
    def summary(self, summaryarg='BB'):
        '''creates .csv files with average results'''
        
        os.makedirs('csv_summary', exist_ok=True)
        
        if summaryarg == 'BB':
            print('Generating csv files with RMSD and reduced ChiSq for backbone atoms ...')
            df = self.ppm_data
            df = df[(df.atomtype.isin(['C','CA','CB','HA','H','N'])) & (df.ChiSq > 0.000)]
            df_grouped_chi = df.groupby(['atomtype'])['ChiSq'].agg({'atomtype':'size', 'ChiSq':'mean'}) \
                .rename(columns={'atomtype':'Entries used','ChiSq':'mean_red_ChiSq'})
            df_grouped_RMSD = (df.groupby(['atomtype'])['error_squared'].mean())**0.5

            df_grouped_chi.to_csv('./csv_summary/output_chi_squared_BB.csv')
            df_grouped_RMSD.to_csv('./csv_summary/output_RMSD_BB.csv')
            
        elif summaryarg == 'CH3_carbons': 
            assert self.ch3_data is not None, 'There is no data from CH3shift'
            print('Generating csv files educed ChiSq for methyl carbons ...')
            df = self.ch3_data
            carbons = ['CD1','CD2','CG2','CD','CG1','CG2','CG2']
            df = df[df.atomtype.isin(carbons) & (df.ChiSq > 0.000)]
            df_grouped_chi = df.groupby(['residue','atomtype'])['ChiSq'].agg({'atomtype':'size', 'ChiSq':'mean'}) \
                .rename(columns={'atomtype':'Entries used','ChiSq':'mean_red_ChiSq'})
            df_grouped_chi.to_csv('./csv_summary/output_chi_squared_CH3.csv')
        
        elif summaryarg == 'ppm_protons':
            print('Generating csv file with protons')
            df = self.ppm_data
            sel = df.atomtype.str.match(r'H[B-Z]+[1-9]*')
            df = df[sel & (df.ChiSq > 0.000)]
            df_grouped_chi = df.groupby(['residue','atomtype'])['ChiSq'].agg({'atomtype':'size', 'ChiSq':'mean'}) \
                .rename(columns={'atomtype':'Entries used','ChiSq':'mean_red_ChiSq'})
            df_grouped_chi.to_csv('./csv_summary/output_chi_squared_H_ppm.csv')
            
                

#    
#
#    
#yy##############################################################################
###############              EXECUTION                        ################
#y##############################################################################
 
    
    

if __name__=='__main__':
#------------------------------------------------------------------------------------------#
#            PIPELINE: reference loading and processing (nomenclature)                     #
#------------------------------------------------------------------------------------------#
        
    ## Choose PDB nomenclature
    pdb_options = {2:'PDBV2',1:'PDBV1'}
    versionPDB = pdb_options[2] #currently only CHARMM27 nomenclature. 
    
        
    if not(bool(args.bmrb_pre) & bool(args.ch3_shift)): #when run just to get the dataframe with experimental data, no need to rename
        reference = args.reference

    else: 
        #transforms nomenclature of the PDB into CHARMM27 nomenclature. Relevant for protons. 
        #No capping groups can be in the reference.
        create_reference_right_nomenclature(args.reference)
        reference = 'reference_renamed.pdb'
        #check nomenclature 
        check_nomenclature(reference_df, versionPDB = versionPDB)

    try:
        referencepdb = PandasPdb().read_pdb(reference)
        reference_df = referencepdb.df['ATOM']
        reference_df = renumber_sequentially(reference_df,'residue_number','residue_name')
        referencepdb.df['ATOM'] = renumber_sequentially(reference_df,'residue_number','residue_name')
 

    except FileNotFoundError as err:
        print(err, '\n', 'The path of the reference PDB is not right. Please, check again. ')
        sys.exit('Execution aborted due to missing reference PDB')
       
    #Pandas dataframe with experimental data
    experimental, seq_exp = experimental_df(args.experimental) #args.experimental
        
    #Renumbering of residues in the experimental dataframe to reference. New numbering goes into column idPDB
    experimental = renumbering_to_reference(reference,experimental) #args.reference

    if not(bool(args.bmrb_pre) & bool(args.ch3_shift)): 
        experimental.to_csv('experimental.csv')



#------------------------------------------------------------------------------------------#
#                                PIPELINE: PPM results                                     #
#------------------------------------------------------------------------------------------#
 

    # checks if there is ppm 
    if bool(args.bmrb_pre):
                
        #Creates pandas data frame for PPM  predictions. 
        predicted_ppm, seq_pre = predicted_ppm_df(args.bmrb_pre) #args.bmrb_pre
        
        #renumbers according to PDB numbering. New numbering goes into column idPDB
        predicted_ppm = renumbering_to_reference(reference,predicted_ppm)
    
    
        #addition of experimental data to 
        left = predicted_ppm
        right = experimental.loc[:,['idPDB','residue','atomtype','CS_exp','error_exp','id_seq']]
        
        
        #df_compiled contains the merged data of experimental and predicted_ppm
        df_merged = pd.merge(left, right, how='outer', on=['idPDB','residue','atomtype'], left_on=None, right_on=None,
             left_index=True, right_index=False, sort=True,
             suffixes=('_x', '_y'), copy=True, indicator=True,
             validate=None)
        df_merged.index = range(len(df_merged))
        df_merged = df_merged[df_merged.idPDB>0].loc[:,['idPDB','residue','atomtype','CS','CS_exp','error_exp','_merge']]
        
        #add_error of PPM 
        df_merged=add_error_ppm(df_merged)
        
        #add extra column with pbd nomenclature
        df_merged['atom_name'] = df_merged.apply(bmrb_to_pdb2,axis=1)
        
        
        #deal with equivalent hydrogens. E.g. HE1, HE2, HE3 appearing as HE in experimental
        #but as HE1,HE2,HE3 in predicted_ppm. This would lead to having HE,HE1,HE2,HE3 
        #in the compiled merged dataframe. We fix it here. 
        df_merged=eliminate_duplicate_hydrogens(df_merged)
        
        #add extra column with pbd nomenclature
        #df_merged['atom_name'] = df_merged.apply(bmrb_to_pdb2,axis=1)
        
        #Error squared and ChiSq
        df_merged['error_squared']=(df_merged['CS']-df_merged['CS_exp'])**2
        df_merged['ChiSq']=df_merged['error_squared']/df_merged['error_predictor']**2
    

#------------------------------------------------------------------------------------------#
#                           PIPELINE: CH3shit results                                      #
#------------------------------------------------------------------------------------------#
 
    
    if bool(args.ch3_shift):#(args.ch3_shift):

        ch3shift=predicted_ch3(args.ch3_shift)
        
        #We discard MET predictions            
        ch3shift=ch3shift[ch3shift.residue!='MET']
        
        #add experimental data
        ch3shift=add_experimental_data_to_ch3shift(ch3shift,experimental)
        
        #Error squared and ChiSq
        ch3shift['error_squared']=(ch3shift['CS']-ch3shift['CS_exp'])**2
        ch3shift['ChiSq']=ch3shift['error_squared']/ch3shift['err2']**2

    else:
        ch3shift=None
    

#------------------------------------------------------------------------------------------#
#                                PIPELINE: PPM output                                      #
#------------------------------------------------------------------------------------------#

    if args.bmrb_pre:
        
        print('Generating PPM output ...') 
        results = visualizer(reference = reference_df, ppm_data = df_merged, ch3_data = None)
    
        print('Generating output in ',os.getcwd())
        # PPM results
        results.summary(summaryarg = 'BB')
        results.summary(summaryarg = 'ppm_protons')
        print(args.pdb_mapping, type(args.pdb_mapping))
        if args.pdb_mapping=='yes':
            results.PDBmapping(plotarg = 'BB')
        if args.histograms=='yes':
            print('doing histograms')
            results.histogramer_plotter()
        df_merged.to_csv('./csv_summary/summary_PPM.csv')

 

#------------------------------------------------------------------------------------------#
#                                PIPELINE: CH3shift output                                 #
#------------------------------------------------------------------------------------------#
 
    if bool(args.ch3_shift):
        print('Generating output from CH3shift ...')

        results = visualizer(reference = reference_df, ppm_data = None, ch3_data = ch3shift)
        results.summary(summaryarg = 'CH3_carbons')

        if args.pdb_mapping=='yes':
            results.PDBmapping(plotarg = 'CH3_carbons')
        ch3shift.to_csv('./csv_summary/summary_CH3shift.csv')





#-------------------------------------------------------------------------------------------
