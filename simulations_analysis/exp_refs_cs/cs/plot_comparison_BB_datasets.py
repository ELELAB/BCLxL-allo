import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 

# trimmed data 
data_18250 = pd.read_csv('summary_PPM_18250.csv')
data_18792 = pd.read_csv('summary_PPM_18792.csv')
data_18793 = pd.read_csv('summary_PPM_18793.csv')


# trimmed data 
exp_CS_18250 = data_18250[['idPDB','residue','atomtype','CS_exp','error_exp']]
exp_CS_18792 = data_18792[['idPDB','residue','atomtype','CS_exp','error_exp']]
exp_CS_18793 = data_18793[['idPDB','residue','atomtype','CS_exp','error_exp']]


# merging datasets 
merged =  pd.merge(exp_CS_18250, exp_CS_18792, how='inner', on=['idPDB','atomtype','residue'], sort=True, suffixes=('_18250', '_18792'),validate=None)
print(merged.columns)
merged = pd.merge(merged, exp_CS_18793, on=['idPDB','atomtype','residue'], sort=True)
print(merged.columns)
merged = merged.rename(columns={'CS_exp':'CS_exp_18793', 'error_exp':'error_exp_18793'})
print(merged.columns)


#generating plots
for atom in ['C','CA','CB','N']: 
    sel = merged[(merged.atomtype == atom) & ((merged.error_exp_18250>0.00) &  (merged.error_exp_18792>0.00) & (merged.error_exp_18793>0.00))]
    plt.figure(figsize=(20,5))
    i#sel = merged[(merged.atomtype == atom) &  ((merged.error_exp_18792>0.00) & (merged.error_exp_18793>0.00))]
    #plots
    plt.errorbar(sel['idPDB'],sel['CS_exp_18250'],sel['error_exp_18250'],fmt='--o',label='18250 (2LPC)')
    plt.errorbar(sel['idPDB'],sel['CS_exp_18792'],sel['error_exp_18792'],fmt='--o',label='18792 (2MO3)')
    plt.errorbar(sel['idPDB'],sel['CS_exp_18793'],sel['error_exp_18793'],fmt='--o',label='18793 (2MO4)')

    plt.legend()
    x_coord = sel['idPDB'].to_list()
    y_coord = np.add(sel['CS_exp_18250'],0.2).to_list()
    
    plt.xlabel('Residue number')
    plt.ylabel('Chemical shift (ppm)')
    plt.title('Atom type {}'.format(atom))
    tag = [str(int(sel['idPDB'].tolist()[j])) + sel['residue'].tolist()[j] for j in range(len(x_coord))] #add tags with residue labels
    
    for i in range(len(x_coord)):    
        plt.text(x_coord[i], y_coord[i], tag[i], fontsize=6)
    plt.savefig('CS_BB_datasets_2_{}.pdf'.format(atom))
    
