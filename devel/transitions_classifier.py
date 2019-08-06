#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 15:34:17 2019

@author: silvia
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from itertools import combinations
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + '/')
from parameters import *
from datetime import datetime





def extract_poly(i1, indexes, poly_cont):
    '''Parses poly.out file extracting all data from each level'''
    if int(i1)+1 == len(indexes):
        starting_line = int(indexes[i1]) + 1
        data = poly_cont[starting_line:-1]
        return data        
    else:
        starting_line = int(indexes[i1]) + 1
        end_line = int(indexes[i1 + 1]) -1
        data = poly_cont[starting_line:end_line]
        return data

def first_d(k1,k2,H):
    '''Calculates first derivative'''
    k1 = float(k1)
    k2 = float(k2)
    H = float(H)
    x = k1 + 2*k2*H
    return x


def classify_crss(H_middle, a, b, index_a, index_b, pend_a, pend_b,AE, Apend, Abs_pend, out_f1, out_f2, thrs_AE, f_comb):
    '''Classifies the type of crossing (avoided, non-avoided) based on pends'''
    f_comb = open(f_comb, 'a')
    thrs_AE2 = float(thrs_AE)/100
    
    if Abs_pend <= thrs_pend:
        if pend_a/pend_b < 0 :
            if pend_a > 0:
                if pend_b < 0:
                    if Apend <= thrs_pend:
                        # Writes
                        H_middle = '{:.7f}'.format(H_middle)
                        a = '{:.7f}'.format(a)
                        b = '{:.7f}'.format(b)
                        pend_a = '{:.7f}'.format(pend_a)
                        pend_b = '{:.7f}'.format(pend_b)
                        AE = '{:.7f}'.format(AE)
                        Apend = '{:.7f}'.format(Apend)
                        Abs_pend = '{:.7f}'.format(Abs_pend)
                
                        f_comb.write(H_middle + '     ' +  a + '     ' + b+ '     ' +AE+ '     '+
                                    pend_a + '     ' + pend_b + '     '  + Abs_pend + '     ' + Apend + '\n')
                        
                        out_f1.write(H_middle + '     ' +  str(index_a) +'     ' +
                                    str(index_b)+ '     ' + a + '     ' + b+ '     ' +
                                    pend_a + '     ' + pend_b + '     ' + '     ' + AE+
                                    '     ' + Apend + '\n')
                
                    elif Apend > 1:
                        if AE <= thrs_AE2:
                            # Writes
                            H_middle = '{:.7f}'.format(H_middle)
                            a = '{:.7f}'.format(a)
                            b = '{:.7f}'.format(b)
                            pend_a = '{:.7f}'.format(pend_a)
                            pend_b = '{:.7f}'.format(pend_b)
                            AE = '{:.7f}'.format(AE)
                            Apend = '{:.7f}'.format(Apend)
                            Abs_pend = '{:.7f}'.format(Abs_pend)
                    
                            out_f2.write(H_middle + '     ' +  str(index_a) +'     ' +
                                        str(index_b)+ '     ' + a + '     ' + b+ '     ' +
                                        pend_a + '     ' + pend_b + '     ' + '     ' + AE+
                                        '     ' + Apend + '\n') 
                            f_comb.write(H_middle + '     ' +  a + '     ' + b+ '     ' +AE+ '     '+
                                    pend_a + '     ' + pend_b + '     '  + Abs_pend + '     ' + Apend + '\n')
            
            elif pend_a < 0:
                if pend_b > 0:
                   if Apend <= thrs_pend:
                        # Writes
                        H_middle = '{:.7f}'.format(H_middle)
                        a = '{:.7f}'.format(a)
                        b = '{:.7f}'.format(b)
                        pend_a = '{:.7f}'.format(pend_a)
                        pend_b = '{:.7f}'.format(pend_b)
                        AE = '{:.7f}'.format(AE)
                        Apend = '{:.7f}'.format(Apend)
                        Abs_pend = '{:.7f}'.format(Abs_pend)
                
                        out_f1.write(H_middle + '     ' +  str(index_a) +'     ' +
                                    str(index_b)+ '     ' + a + '     ' + b+ '     ' +
                                    pend_a + '     ' + pend_b + '     ' + '     ' + AE+
                                    '     ' + Apend + '\n')
                        f_comb.write(H_middle + '     ' +  a + '     ' + b+ '     ' +AE+ '     '+
                                    pend_a + '     ' + pend_b + '     '  + Abs_pend + '     ' + Apend + '\n')
                
                   elif Apend > 1:
                         if AE <= thrs_AE2:
                            # Writes
                            H_middle = '{:.7f}'.format(H_middle)
                            a = '{:.7f}'.format(a)
                            b = '{:.7f}'.format(b)
                            pend_a = '{:.7f}'.format(pend_a)
                            pend_b = '{:.7f}'.format(pend_b)
                            AE = '{:.7f}'.format(AE)
                            Apend = '{:.7f}'.format(Apend)
                            Abs_pend = '{:.7f}'.format(Abs_pend)
                    
                            out_f2.write(H_middle + '     ' +  str(index_a) +'     ' +
                                        str(index_b)+ '     ' + a + '     ' + b+ '     ' +
                                        pend_a + '     ' + pend_b + '     ' + '     ' + AE+
                                        '     ' + Apend + '\n') 
                            f_comb.write(H_middle + '     ' +  a + '     ' + b+ '     ' +AE+ '     '+
                                    pend_a + '     ' + pend_b + '     '  + Abs_pend + '     ' + Apend + '\n')
        elif pend_a/pend_b > 0:
            print 'exotic crosses'
    f_comb.close()


def limits_avoided(k0_a, k0_b, k1_a, k1_b, k2_a, k2_b):
    '''Calculates first derivative of an avoided crossing'''
    a = float(k2_a)
    b = float(k1_a)
    c = float(k0_a)
    a2 = float(k2_b)
    b2 = float(k1_b)
    c2 = float(k0_b)
    
    x1 = -(b)/(2*(a))
    x2 = -(b2)/(2*(a2))
    y1 = c - ((b**2)/(4*a))
    y2 = c2 - ((b2**2)/(4*a2))
    return x1, x2, y1, y2


def limits_nonavoided(k0_a, k0_b, k1_a, k1_b, k2_a, k2_b):
    '''Calculates crossing H value of a non-avoided crossing'''
    a = float(k2_a - k2_b)
    b = float(k1_a - k1_b)
    c = float(k0_a - k0_b)
    sol = []
     
    if a != 0:
        x1 = (-b + math.sqrt(b**2 - 4*a*c)) / (2 * a)
        x2 = (-b - math.sqrt(b**2 - 4*a*c)) / (2 * a)
        print 'Soluciones de la ecuacion: x1=%4.3f y x2=%4.3f ' % (x1, x2)
        sol.append(x1)
        sol.append(x2)
    else:
        if b != 0:
           x = -c / b
           print 'Solucion de la ecuacion: x=%4.3f ' % x
           sol.append(x)
        else:
           if c != 0:
              print 'La ecuacion no tiene solucion. '
           else:
              print 'La ecuacion tiene infinitas soluciones. '
    return sol


def calc_curvature(xmax, k1_a, k1_b, k2_a, k2_b):
    xmax, b1, b2, a1, a2 = float(xmax), float(k1_a), float(k1_b), float(k2_a), float(k2_b)
    curv1 = abs(2*a1)/(1+((2*a1*xmax) + b1)**2)**float(3)/2
    curv2 = abs(2*a2)/(1+((2*a2*xmax) + b2)**2)**float(3)/2
    print 'curv1: ' + str(curv1)
    print 'curv2: ' + str(curv2)
    return curv1, curv2



if __name__ == "__main__":

    
    #
    ##
    ### Reads and parses poli.txt file
    ##
    #
    # Opens and reads
    in_f = open(path + 'res/poli.out', 'r')
    poly_cont = in_f.read().split('\n')
    in_f.close()
    
    # Extracts number of the starting line of each level in ene file 
    indexes = []
    c0 = 0
    for line in poly_cont:
        c0 +=1
        if 'lvl_' in line:
            indexes.append(int(c0)-1)
    in_f.close()    
    
    # Extract poli data
    dic = { i: { H:[] for H in H_values } for i in lvls_list }
    
    # Fill dic with data
    for i1 in range(len(indexes)):
        index_lvl = int(i1+1)
        name_lvl = 'lvl_' + str(index_lvl)
        cont = extract_poly(i1, indexes, poly_cont)
        
        for i2 in range(len(cont)):
            H_value = H_values[i2+4]
            vec = cont[i2].split()
            if len(vec) > 1:
                H = float(vec[1])
                k0 = float(vec[4])
                k1 = float(vec[3])
                k2 = float(vec[2])
                dic[name_lvl][H] = [k0, k1, k2] 
            else:
                continue
            



    # Load input file
    final_df = pd.read_csv(path + 'res/final.txt', sep='\t')
    
    # Open output files
    out_avoid = open('res/avoided_sum.txt', 'w')
    out_nonavoid = open('res/nonavoided_sum.txt', 'w')
    
    
    # Check if subdir_path exists and if its empty
    subdir_path = path + 'combinations/'
    if not os.path.exists(subdir_path):
        os.makedirs(subdir_path)
    
    else:
        ldir = os.listdir(subdir_path)
        if len(ldir) > 0:
            for i in ldir:
                os.remove(subdir_path + i)
    
        
    #
    ##
    ### Main loop
    ##
    #
    
    H_s = H_values.astype(str)
    dic_2 = { H: { i:[] for i in lvls_list } for H in H_values }
    
    for i in range(dim[0]):
        
        if i+4 >= dim[0]-1:
            break    
        
        else:
            H_value = H_values[i+4]
            H_middle = H_values[i+2]
            sEs = final_df.iloc[i+2]
            print final_df
            
            for a,b in combinations(sEs,2):  ### campo  middle
                # Calculates AE
                AE = abs(a-b)
                
                #print sEs
                #Retrieves index level
                index_a = sEs[sEs == a].index[0]
                index_b = sEs[sEs == b].index[0]
                
                
                # Output file
                f_comb = subdir_path + 'comb_' + index_a + '_' + index_b + '.out'
                print f_comb
                
                # Checks if file already exists
                if not os.path.isfile(f_comb):
                    # Creates file
                    f = open(subdir_path + 'comb_' + index_a + '_' + index_b + '.out' , 'w+')
                    # Writes data in file
                    f.write(u'#H   Ei   Ej   AEij   pendi   pendj   |pendi|-|pendj|   |pendi-pendj|\n')
                    f.close()
                
                
                if AE <= thrs_AE:
                    
                    # Saves levels somewhere, to then evaluate change in signs
                    k1_a = float(dic[index_a][H_value][1])
                    k1_b = float(dic[index_b][H_value][1])
                    k2_a = float(dic[index_a][H_value][2])
                    k2_b = float(dic[index_b][H_value][2])
                    
                    # Calculates pend
                    pend_a = first_d(k1_a, k2_a, H_middle)
                    pend_b = first_d(k1_b, k2_b, H_middle)
                    
                    # Pend differences
                    Apend = abs(pend_a - pend_b)
                    Abs_pend = abs(pend_a) - abs(pend_b) 
                    div_pend = pend_a / pend_b
                    
                    # Classify non-avoided // avoided crossings
                    classify_crss(H_middle, a, b, index_a, index_b, pend_a, pend_b, AE, Apend, Abs_pend, out_avoid, out_nonavoid, thrs_AE, f_comb)
            
    out_avoid.close()
    out_nonavoid.close()
    
    
    #
    ## Post-analysis of crosses
    #
    ldir = os.listdir(subdir_path)
    res = open(path + 'res/results.out', 'w')
    
    # Write header res file
    res.write('#Level_a     Level_b     Type_crossing     H_crossing     E_crossing1     E_crossing2     AE     Curvature1     Curvature2\n')
    
    # Loop
    for f in ldir:
        if not os.stat(subdir_path + f).st_size > 71:
            os.remove(subdir_path+f)
        
        else:
            # Read file
            name = f.strip('comb_').strip('.out').split('_',2)
            a = name[0] +'_' +  name[1]
            b = name[2]
            array = np.loadtxt(subdir_path+f, comments = '#')
            dim = np.shape(array)
            H = array[:,0]
            E_a = array[:,1] 
            E_b = array[:,2]
            AEab = array[:,3]
            pend_a = array[:,4]
            pend_b = array[:,5]
            diff1 = array[:,6]
            diff2 = array[:,7]
            
            sign_penda = np.sign(array[:,4])
            sign_pendb = np.sign(array[:,5])
            
            signa = np.where(np.diff(np.sign(array[:,4])))[0]
            signb = np.where(np.diff(np.sign(array[:,5])))[0]
            
            ordered_AE = np.sort(array[:,3])
            minimAE = ordered_AE[0]
            minimAE_2 = ordered_AE[1]
            i1, j1 = np.where(array == minimAE)
            i2, j2 = np.where(array == minimAE_2)
            H_min = array[i1, 0]
            H_min2 = array[i2, 0]
            
            # Find polynomio for Hmin + 2 
            H_values = np.array(H_values)
            H_poly_ind = np.where(H_values == H_min[0])[0][0]
            H_poly_ind2 = np.where(H_values == H_min2[0])[0][0]
            
            H_poly = H_values[H_poly_ind + 2]
            H_poly2 = H_values[H_poly_ind2 + 2]
            
            # Assign data
            k0_a = dic[a][H_poly][0]
            k1_a = dic[a][H_poly][1]
            k2_a = dic[a][H_poly][2]
            k0_b = dic[b][H_poly][0]
            k1_b = dic[b][H_poly][1]
            k2_b = dic[b][H_poly][2]
            
            
            if len(signa) < 1:
                if len(signb) <1:
                    # Non-avoided crossings
                    type_cross = 'non-avoided'
                    print '1: ' + a + ' ' + b
                    # Calculates limits
                    res_v =  limits_nonavoided(k0_a, k0_b, k1_a, k1_b, k2_a, k2_b)
                    
                    # Select the good solution
                    if len(res_v) > 1:
                        if abs(H_min - res_v[0]) < abs(H_min - res_v[1]):
                            H_cross = res_v[0]
                        elif abs(H_min - res_v[0]) > abs(H_min - res_v[1]):
                            H_cross = res_v[1]
                    else:
                        H_cross = res_v[0]
                    
                    # Calculates E_cross
                    E_cross = (k2_a * (H_cross**2)) + (k1_a * H_cross) + k0_a
                    # Writes output
                    H_cross = '{:.5f}'.format(H_cross)
                    E_cross = '{:.5f}'.format(E_cross)
                    res.write(a + '     '+ b + '     '+ type_cross + '     '+ str(H_cross) + '     '+ str(E_cross) + '     -     - \n' )
                    
                    
                else:
                    # Avoided crossing
                    type_cross = 'avoided'
                    print '2: '+ a + b
                    res_v = limits_avoided(k0_a, k0_b, k1_a, k1_b, k2_a, k2_b)
                    H_cross = res_v[0]
                    E_cross1 = res_v[2]
                    E_cross2 = res_v[3]
                    AE = abs(E_cross1 - E_cross2)
                    AE = '{:.5f}'.format(AE)
                    H_cross = '{:.5f}'.format(H_cross)
                    E_cross1 = '{:.5f}'.format(E_cross1)
                    E_cross2 = '{:.5f}'.format(E_cross2)
                    curvature = calc_curvature(H_cross, k1_a, k1_b, k2_a, k2_b)
                    curv1 = '{:.5f}'.format(curvature[0])
                    curv2 = '{:.5f}'.format(curvature[1])
                    #curvature = '{:.5f}'.format(curvature)
                    res.write(a + '     '+ b + '     '+ type_cross + '     '+ H_cross + '     '+
                              E_cross1 + '     ' + E_cross2 + '     ' + AE +'     ' + curv1 +'     ' + curv2 +  '\n' )
                    
                    
            else: # a es linea curva
                if len(signb) < 1:
                    # Non-avoided crossings
                    print '3: '+ a + b
                    res_v = limits_nonavoided(k0_a, k0_b, k1_a, k1_b, k2_a, k2_b)
                    # Calculates limits
                    res_v =  limits_nonavoided(k0_a, k0_b, k1_a, k1_b, k2_a, k2_b)
                    # Check results
                    if len(res_v) > 1:
                        if abs(H_min - res_v[0]) < abs(H_min - res_v[1]):
                            H_cross = res_v[0]
                        elif abs(H_min - res_v[0]) > abs(H_min - res_v[1]):
                            H_cross = res_v[1]
                    else:
                        H_cross = res_v[0]
                    # Calculates E_cross
                    E_cross = (k2_a * (H_cross**2)) + (k1_a * H_cross) + k0_a
                    H_cross = '{:.5f}'.format(H_cross)
                    E_cross = '{:.5f}'.format(E_cross)
                    # Writes output
                    res.write(a + '     '+ b + '     '+ type_cross + '     '+ str(H_cross) + '     '+ str(E_cross) + '     -     - \n' )
                    
                else:
                    # Avoided crossing
                    type_cross = 'avoided'
                    print '4: '+ a + b                
                    res_v = limits_avoided(k0_a, k0_b, k1_a, k1_b, k2_a, k2_b)
                    H_cross = res_v[0]
                    E_cross1 = res_v[2]
                    E_cross2 = res_v[3]
                    AE = abs(E_cross1 - E_cross2)
                    AE = '{:.5f}'.format(AE)
                    H_cross = '{:.5f}'.format(H_cross)
                    E_cross1 = '{:.5f}'.format(E_cross1)
                    E_cross2 = '{:.5f}'.format(E_cross2)
                    curvature = calc_curvature(H_cross, k1_a, k1_b, k2_a, k2_b)
                    curv1 = '{:.5f}'.format(curvature[0])
                    curv2 = '{:.5f}'.format(curvature[1])
                    res.write(a + '     '+ b + '     '+ type_cross + '     '+ H_cross + '     '+
                              E_cross1 + '     ' + E_cross2 + '     ' + AE + '     ' + curv1 +'     ' + curv2 +  '\n' )
                    
    
    res.close()
        
    
    
    

    print "#################################################"
    print 'End of run'
    print "#################################################"
    print ''
    print "Check you results at /res/ directory"
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print 'Bye!'