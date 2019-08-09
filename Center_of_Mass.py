#!/usr/bin/env python
# coding: utf-8

# In[1]:


#Description: Calculate center of mass of a protein from it's pdb ID
#Author: Rik Dhar


# In[2]:


import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
ppdb = PandasPdb()

import warnings
warnings.filterwarnings("ignore")

#Importing library pandas, numpy and biopandas


# In[14]:


name = input("Enter PDB ID: ")
type(name)
ppdb.fetch_pdb(name)

#Asking user input for pdb ID and then fetching file as a pandas dataframe.


# In[4]:


hetatm_df = ppdb.df['HETATM']
atm_df = ppdb.df['ATOM']
het_woH2O_df = hetatm_df[hetatm_df["residue_name"].str.contains('HOH') == False]

#Creating a new dataframes for atoms, hetroatoms, and heteroatoms without water.


# In[5]:


nan = np.nan
atomic_mass = pd.DataFrame({'symbol': ['H','HE','LI','BE','B','C','N','O','F','NE','NA','MG','AL','SI','P','S','CL','K','AR','CA','SC','TI','V','CR','MN','FE','NI','CO','CU','ZN','GA','GE','AS','SE','BR','KR','RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG','CD','IN','SN','SB','I','TE','XE','CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY','HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT','AU','HG','TL','PB','BI','PO','AT','RN','FR','RA','AC','PA','TH','NP','U','AM','PU','CM','BK','CF','ES','FM','MD','NO','RF','LR','DB','BH','SG','MT','RG','HS','DS','CN','NH','FL','MC','LV','TS','OG'],
                            'atomic mass': [1.0079,4.0026,6.941,9.0122,10.811,12.0107,14.0067,15.9994,18.9984,20.1797,22.9897,24.305,26.9815,28.0855,30.9738,32.065,35.453,39.0983,39.948,40.078,44.9559,47.867,50.9415,51.9961,54.938,55.845,58.6934,58.9332,63.546,65.39,69.723,72.64,74.9216,78.96,79.904,83.8,85.4678,87.62,88.9059,91.224,92.9064,95.94,98,101.07,102.9055,106.42,107.8682,112.411,114.818,118.71,121.76,126.9045,127.6,131.293,132.9055,137.327,138.9055,140.116,140.9077,144.24,145,150.36,151.964,157.25,158.9253,162.5,164.9303,167.259,168.9342,173.04,174.967,178.49,180.9479,183.84,186.207,190.23,192.217,195.078,196.9665,200.59,204.3833,207.2,208.9804,209,210,222,223,226,227,231.0359,232.0381,237,238.0289,243,244,247,247,251,252,257,258,259,261,262,262,264,266,268,272,277,nan,nan,nan,nan,nan,nan,nan,nan,]})  

#Creating a dataframe with the list of atomic masses


# In[6]:


def mass_multiply(df):
    df['mass'] = nan
    for x in range(df.index.size):
        X = df.at[x,'element_symbol']
        for y in range(atomic_mass.index.size):
            Y = atomic_mass.at[y,'symbol']
            if X == Y:
                df.at[x,'mass'] = atomic_mass.at[y,'atomic mass'] 
    df['x_coord*mass'] = df.x_coord*df.mass
    df['y_coord*mass'] = df.y_coord*df.mass
    df['z_coord*mass'] = df.z_coord*df.mass 

#Function to create new columns in the dataframe with the coordinates multiplied with their respective masses.


# In[7]:


mass_multiply(atm_df)
mass_multiply(hetatm_df)
mass_multiply(het_woH2O_df)

#Applying function to different dataframes


# In[8]:


def com(df):
    x = df['x_coord*mass'].sum()/df['mass'].sum() 
    y = df['y_coord*mass'].sum()/df['mass'].sum()
    z = df['z_coord*mass'].sum()/df['mass'].sum()
    return {'CoM_X': x, 'CoM_Y': y ,'CoM_Z': z}

#Function to calculate center of mass and return the coordinates as a dictionary


# In[9]:


merge = [atm_df, het_woH2O_df]
total_df = pd.concat(merge)

#Merge the atom dataframe with the hetero atom with water dataframe


# In[10]:


print("Center of mass of atoms of the protein:")
print(com(atm_df))
print("\n")
print("Center of mass of heteroatoms of the protein:")
print(com(hetatm_df))
print("\n")
print("Center of mass of heteroatoms of the protein without water:")
print(com(het_woH2O_df))
print("\n")
print("Center of mass of atoms and the heteroatoms without water:")
print(com(total_df))


# hetatm_df.to_csv(r'C:\Users\rikdh\Documents\Lab Work\Computational\Project_Center_of_Mass\HETATOM.csv')

# string = "hetatm_df['x_coord*mass'] = hetatm_df.x_coord*hetatm_df.masshetatm_df['y_coord*mass'] = hetatm_df.y_coord*hetatm_df.masshetatm_df['z_coord*mass'] = hetatm_df.z_coord*hetatm_df.mass"
# string = string.replace('hetatm_df', 'type')
# hetatm_df.loc[0] 
# string

# In[ ]:




