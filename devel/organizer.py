#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 12:29:04 2019

@author: silvia
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/')
from parameters import *
from datetime import datetime
from matplotlib.gridspec import GridSpec





def poly(five_H, five_E): #dataset with set of points: would retrieve 
    '''Calculates the polynomia expression for a set of H and E'''
    x = five_H # set campos H (5)
    y = five_E # set valores E (5)
    vec_z = np.polyfit(x, y, deg = 2)
    return vec_z


def expected_E(vec_z, H_v):
    '''Retrieves the E value of a spin level at a given H'''
    k0 = vec_z[2]
    k1 = vec_z[1]
    k2 = vec_z[0]
    predicted_E = (k2 * (float(H_v)**2)) + (k1 * float(H_v)) + k0
    return predicted_E
    


def search(v0, v1, v2, list_e): #index = pos-1
    '''Finds the corresponding E value in simpre.ene comparing with expected_value'''
    #num = []
    diffs = {}
    for i in range(len(list_e)):
        e0 = list_e[i]
        if e0 == 'u':
            continue
        
        
        e0 = round(float(e0), 7)
        diff = abs(float(v0) - float(e0))
        diffs[i] = diff
        
        """
        if  diff < 10**-6:
            num.append(i)
        else:
            if  e0 <= v1:
                if e0 >= v2:
                    num.append(i)
        """
    num = min(diffs, key=diffs.get)
    return num







if __name__ == "__main__":    
    #
    ##
    ### Main loop
    ##
    #
    # Open out files
    poli_out = open('res/poli.out', 'w')
    
    col = 0
    row = 0
    for i0 in range(dim[1]):
            
        five_H = []
        five_E = []
        d_poli = {}
        
        poli_out.write('{0:<2}{1}\n'.format('#', lvls_list[i0]))
        poli_out.write('{0} {1:^10} {2:^10} {3:^10} {4:^10}\n'.format('#', 'H' , 'k2', 'k1', 'k0'))
                       
        for i1 in range(dim[0]):
            
            if i1 == dim[0]-1:
                print i1
                print dim[0]
                continue
            
            # Valores de campo : actual y next
            H_value = H_values[i1]
            H_next = H_values[i1+1]
            
            # Rellenar lista H
            five_H.append(H_value)
            
            # Rellenar lista E
            #numero de columna de ene_df correspondiente al nivel i0 para el campo H i1
            pos = order_df.iloc[i1, i0]
            if pos == 'nan':
                five_E.append(float(expected_df.iloc[i1, i0]))
            
            else:
                pos = int(pos)
                five_E.append(ene_df.iloc[i1, pos-1])
            
            # Everything starts now......
            if len(five_H) ==5:
                            
                # Fitting
                vec_poly = poly(five_H, five_E)
                p = np.poly1d(vec_poly)
                d_poli[i1] = p
                vec = ['{:.7f}'.format(i) for i in vec_poly]
                poli_out.write('{0:>10} {1:>10} {2:>10} {3:>10}\n'.format(str(H_value), vec[0], vec[1], vec[2]))
                
                # Calculates expected_v
                v0 = str(round(float(expected_E(vec_poly, H_next)), 7))
                v1 = float(v0) + half_minAE
                v2 = float(v0) - half_minAE
                v1 = round(v1,7)
                v2 = round(v2,7)
                
                # Stores expected value of level i0 at next_H value
                expected_df.iloc[i1+1, i0] = v0
                
                # Searches in list of Es from ene at one H value (i1+1)
                s_es = (ene_df.iloc[i1+1]).tolist()
                num = search(v0, v1, v2, s_es)
                order_df.iloc[i1+1, i0] = int(num)+1
                
                
                five_E.pop(0)
                five_H.pop(0)
            
        poli_out.write('\n')
    
    # Write outputs
    order_df.to_csv(path + 'res/order.txt', header = lvls_list, sep='\t', na_rep='na')
    expected_df.to_csv(path + 'res/expected.txt', header = lvls_list, sep='\t', na_rep='na')
    poli_out.close()
    
    #
    ##
    ### Obtain final_df
    ##
    #
    newene_df = np.array(ene_df)
    order_df = np.array(order_df)
    
    dims = ene_df.shape
    tsize = dims[0]*dims[1]
    
    #Creates the final_df array using the labels from order_df as a mask over ene_df
    df_final = np.arange(tsize).reshape(dims[0], dims[1])
    df_final = np.array(df_final, dtype='float32')
    
    for i in range(dims[0]): # valores de H
        for j in range(dims[1]): # niveles
            i = int(i)
            j = int(j)
            
            if order_df[i,j] == 'nan':
                final_df.iloc[i,j] = newene_df[i,j]
            else:
                col = int(order_df[i,j])
                final_df.iloc[i,j] = newene_df[i,col-1] 
            
    
    # Writing outputs
    final_df.to_csv(path + 'res/final.txt', header = lvls_list, sep='\t', na_rep='na', float_format='%.7f' )



# =============================================================================
#   Quality control of the results
# =============================================================================
    # Difference between expected vs. final order
    t = np.empty([dims[0], dims[1]], dtype=float)
    for i0 in range(dims[1]):
        for i1 in range(dims[0]):
            v = abs(float(expected_df.iloc[i1, i0]) - float(final_df.iloc[i1, i0]))
            t[i1,i0] = v    
    
    # Check if subdir_path exists and if its empty
    plotsdir_path = path + 'plots/'
    if not os.path.exists(plotsdir_path):
        os.makedirs(plotsdir_path)
    
    # Visual checking of results final vs. expected
    fig = plt.figure(constrained_layout = True, figsize=(7, 5))
    #fig.suptitle("Results summary")
    fig.text(-0.05, 0.5, 'Energy (cm$^{-1}$)', va='center', rotation='vertical', fontsize = 13 )
    fig.text(0.5,-0.05, 'Magnetic field (T)', ha='center', fontsize = 13)
    
    gs = GridSpec(3, 3, figure=fig, height_ratios=[5, 1,1], width_ratios=[1, 1,1]) #, left=0.05, right=0.48, wspace=0.05)
    
    ax1 = fig.add_subplot(gs[0,0])
    ax1.plot(ene[:,0], ene_df.iloc[:,0:])
    ax1.set_title('Raw data, unsorted')
    
    ax2 = fig.add_subplot(gs[0,1])
    ax2.plot(expected_df.index, expected_df.iloc[:,0:])
    ax2.set_title('Polynomial fit ')
    ax2.tick_params(labelbottom=True, labelleft=False)
    
    ax3 = fig.add_subplot(gs[0,2])
    ax3.plot(final_df.index, final_df.iloc[:,0:])
    ax3.set_title('Raw data, sorted')
    ax3.tick_params(labelbottom=True, labelleft=False)
    
    ax4 = fig.add_subplot(gs[1, :])
    ax4.plot(expected_df.index, t[:,:])
    ax4.set_ylim([-0.25,0.25])
    ax4.set_title('Quality control of expected energies')
    ax4.tick_params(labelbottom=False, labelleft=True)
    
    ax5 = fig.add_subplot(gs[2, :])
    ax5.plot(expected_df.index, t[:,:])
    ax5.set_yscale('log')

    plt.savefig('plots/summary_plots.png', dpi = 300, bbox_inches='tight')
    


# =============================================================================
#   End of run
# =============================================================================
    
    print "#################################################"
    print 'End of run'
    print "#################################################"
    print "Check the quality control at /plots/summary_plots.png"
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print 'Bye!'
    