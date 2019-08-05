#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 12:32:50 2019

@author: silvia
"""

import os
import numpy as np
import pandas as pd


    
# =============================================================================
# DEFINITIONS
# Define here the relevant parameters of your system
# =============================================================================

# Working path
path = os.getcwd() + '/'
#print path

script_path = os.path.realpath(__file__)
#print __file__
#print os.path.join(os.path.dirname(__file__), '..')
#print os.path.dirname(os.path.realpath(__file__))
#print os.path.abspath(os.path.dirname(__file__))

# Check if res_path exists otherwise creates it
resdir_path = path + 'res/'
if not os.path.exists(resdir_path):
    os.makedirs(resdir_path)



# Electronic and nuclear spin
J = float(8)
I = float(7.5)
#J = float(raw_input('Valor de J: '))
#I = float(raw_input('Valor de I: '))
g_par = 0.0000
g_per = 0.0000

#A_par = 0
#A_per = 0

#n_levels = raw_input('How many n lower-energy levels are you interested in?:')
n_levels = 16

# Total spin levels of the system
tot_lvls = int((2*J + 1)*(2*I + 1))
m_j = list(range(int(-(2*J)), int((2*J+2)), int(2)))
m_i = list(range(int(-(2*I)), int((2*I+2)), int(2)))





if __name__ == "__main__":

    ## Generate projections list from m_j and m_i
    projections_l = []
    for j in m_j:
        for i in m_i:
            proj = str(j) + '/' + str(i)
            projections_l.append(proj)
    
    # Loads new 'simpre.ene' file with energies of each level (quantum number) at diff H field
    with open(path + 'simpre.ene', 'r') as ene_f:
        ene = np.loadtxt(ene_f, dtype= float)
        if n_levels == '':
            n_levels = np.size(ene, axis=1)-1
            ene_df = pd.DataFrame(data=ene[0:,1:],  index=ene[0:,0]) #, columns=range(1,len(lvls_list) +1))
        else:
            n_levels = int(n_levels)
            ene_df = pd.DataFrame(data=ene[0:,1:(n_levels + 1)],  index=ene[0:,0])
    
    
    # Removes values at 0.0T
    if float(0.0) in ene_df.index:
        ene_df = ene_df.drop([0.0])
    
    
    
    # Generates levels list    
    lvls_list = []
    for i in list(range(n_levels)):
        lvls_list.append('lvl_'+ str(i+1))
    
    
    
    
    
    
       
    # Creating an empty DF for further lvl ordering
    H_values = ene_df.index
    expected_df = pd.DataFrame(index=H_values, columns=lvls_list)
    final_df = pd.DataFrame(index=H_values, columns=lvls_list)
    order_df = pd.DataFrame(index=H_values, columns=lvls_list)
    
    # Fill first 5 positions of order_df
    c = 0
    for i in lvls_list:
        c=0
        for i2 in H_values:
            c+=1
            if c < 6 :
                order_df.loc[i2, i] = i.replace('lvl_' , '')
    
    
    
    ## First check
    dim = ene_df.shape
    if int(dim[1]) !=  len(lvls_list):
        print 'hay un problema de dimensiones'
    
    
    
    # Get minimum value 
    minAE = 10000000000
    for i1 in range(dim[0]-1):
        for i2 in range(dim[1]-1):
            v1 = ene_df.iloc[i1, i2]
            v2 = ene_df.iloc[i1, i2+1]
            diff_0 = abs(v1 - v2)
            if diff_0 < minAE:
                minAE = diff_0
    half_minAE = minAE/2