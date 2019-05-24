#!/usr/bin/env python
# coding: utf-8

# Imports and path definitions

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from itertools import combinations

path = os.getcwd() + '/'


# # NEW APPROACH

# #### Pre analysis step!
# 
# -- hiperfino = 0 
# -- Chequear el buen control sobre niveles
# 
# #### First trials:
# ~ Parallel and isotrope hyperfine
# Polinomios de intervalos  5/3 puntos de H
# 
# 
# # Flow
# 1. Run your simpre calculation and plot your results.
# 2. Define a trainning data set region without non-avoided crossings
#     - Minimum of 3-5 points (plot your results)
# 3. Run Poly.py 
# 4. Run expected_values.py
# 5. Run main.py

# In[2]:


ene_f = path + "simpre.ene"

ene = np.loadtxt(ene_f, dtype= float)

plt.plot(ene[:,0], ene[:,[1,2,3,4,5,6,7,8]])
#plt.plot(ene[:,0], ene[:,[17, 18,19,20, 21, 22, 23, 24,25,26]])#, 9,10, 11, 12, 13, 14, 15, 16]])
#plt.plot(ene[:,0], ene[:,1:])
#plt.savefig('plot_ene.png', dpi = 300)

print 'n. of points: ' + str(len(ene))
plt.show()


# In[14]:

# Evaluate expected values

plt.plot(expected_df.index, expected_df.iloc[:,1:])
#plt.plot(expected_df.index, expected_df.iloc[:,[8,9,10,]])
#plt.plot(expected_df.index, expected_df.iloc[:,[16,17, 18,19,20,21,22,23,24,25]])
#plt.savefig('plot_prev.png', dpi = 300)

#print 'n. of points: ' + str(len(ene_prev))
plt.show()
print expected_df

# In[15]

#
## Functions module
#

def poly(five_H, five_E): #dataset with set of points: would retrieve 
    x = five_H # set campos H (5)
    y = five_E # set valores E (5)
    vec_z = np.polyfit(x, y, deg = 2)
    return vec_z



def expected_E(vec_z, H_v):
    k0 = vec_z[2]
    k1 = vec_z[1]
    k2 = vec_z[0]
    predicted_E = (k2 * (float(H_v)**2)) + (k1 * float(H_v)) + k0
    return predicted_E
    


def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])




def search(v0, v1, v2, list_e): #index = pos-1
    num = []
    
    for i in range(len(list_e)):
        e0 = list_e[i]
        if e0 == 'u':
            continue
        
        e0 = round(float(e0), 7)
        diff = abs(abs(float(v0)) - abs(float(e0)))
        if  diff < 10**-6:
            num.append(i)
        else:
            if  e0 <= v1:
                if e0 >= v2:
                    num.append(i)
    return num


# In[4]:


#
## Initial definitions
#

# Defined by the user so far; automatize for the future
J = float(8)
I = float(7.5)
#J = float(raw_input('Valor de J: '))
#I = float(raw_input('Valor de I: '))
g_par = 0.0000
g_per = 0.0000


# Lists
tot_lvls = int((2*J + 1)*(2*I + 1))
m_j = list(range(int(-(2*J)), int((2*J+2)), int(2)))
m_i = list(range(int(-(2*I)), int((2*I+2)), int(2)))

lvls_list = []
for i in list(range(np.size(ene, axis=1)-1)):
    lvls_list.append('lvl_'+ str(i+1))


## Generate projections list from m_j and m_i
projections_l = []
for j in m_j:
    for i in m_i:
        proj = str(j) + '/' + str(i)
        projections_l.append(proj)



# Loads new 'simpre.ene' file with energies of each level (quantum number) at diff H field
with open(path + 'simpre.ene', 'r') as ene_f:
    ene = np.loadtxt(ene_f, dtype= float)
    ene_df = pd.DataFrame(data=ene[0:,1:],  index=ene[0:,0]) #, columns=range(1,len(lvls_list) +1))

if float(0.0) in ene_df.index:
    ene_df = ene_df.drop([0.0])
    


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


# In[]:
#
##
### Main loop
##
#

# Open out files
poli_out = open('poli.out', 'w')

col = 0
row = 0


for i0 in range(dim[1]):
        
    five_H = []
    five_E = []
    d_poli = {}
    
    poli_out.write(lvls_list[i0] + '\n' )
    
    
    for i1 in range(dim[0]):
        
        if i1 == dim[0]-1:
            print i1
            print dim[0]
            continue
        
        
        # Valores de campo : actual y next
        H_value = H_values[i1]
        H_next = H_values[i1+1]
        
        print 'valor de campo actual'
        print H_value
        print 'valor de campo next'
        print H_next
        
        # Rellenar lista H
        five_H.append(H_value)
        
        # Rellenar lista E
        #numero de columna de ene_df correspondiente al nivel i0 para el campo H i1
        pos = int(order_df.iloc[i1, i0])
        five_E.append(ene_df.iloc[i1, pos-1])
       
        
        
        # Everything starts now......
        if len(five_H) ==5:
            
            print 'aqui los 5 campos'
            print five_H
            print 'aqui las 5 Es'
            print five_E
            
            
            # Fitting
            vec_poly = poly(five_H, five_E)
            p = np.poly1d(vec_poly)
            d_poli[i1] = p
            poli_out.write('H: ' + str(H_value) + ' ' + str(vec_poly) +'\n') # Last H value of the poly
            
            
            # Calculates expected_v
            v0 = str(round(float(expected_E(vec_poly, H_next)), 7))
            print 'aqui expected value'
            print v0
            v1 = float(v0) + half_minAE
            v2 = float(v0) - half_minAE
            v1 = round(v1,7)
            v2 = round(v2,7)
            
            
            # Stores expected value of level i0 at next_H value
            expected_df.iloc[i1+1, i0] = v0
            
            # Search in list of Es from ene at one H value (i1+1)
            s_es = (ene_df.iloc[i1+1]).tolist()
            
            num = search(v0, v1, v2, s_es)
            print num
            
            
            if len(num) == 1:
                order_df.iloc[i1+1, i0] = int(num[0])+1
            else:
                print 'ninguno coincide'
                pass
                
        #       print 'isnull'
        #       five_E.append(v0)
        #       print five_E


                

            five_E.pop(0)
            five_H.pop(0)
    
    poli_out.write('\n')

order_df.to_csv(path + 'order.txt', header = lvls_list, sep=' ', na_rep='na')
expected_df.to_csv(path + 'expected.txt', header = lvls_list, sep=' ', na_rep='na')
poli_out.close()
    

# In[]:
