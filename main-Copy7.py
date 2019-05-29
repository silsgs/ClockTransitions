#!/usr/bin/env python
# coding: utf-8

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

# In[1]:


# Imports and path definitions

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from itertools import combinations

path = os.getcwd() + '/'


# In[12]:


ene_f = path + "simpre.ene"

ene = np.loadtxt(ene_f, dtype= float)

#plt.plot(ene[:,0], ene[:,[1,2,3,4,5,6,7,8, 9,10,11,12,13,14]])
#plt.plot(ene[:,0], ene[:,[17, 18,19,20, 21, 22, 23, 24,25,26]])#, 9,10, 11, 12, 13, 14, 15, 16]])
plt.plot(ene[:,0], ene[:,1:])

#plt.savefig('plot_ene.png', dpi = 300)

print 'n. of points: ' + str(len(ene))
plt.show()

# In[12]:


final_f = path + "final.txt"

fin = np.loadtxt(final_f, dtype= float)

#plt.plot(ene[:,0], ene[:,[1,2,3,4,5,6,7,8, 9,10,11,12,13,14]])
#plt.plot(ene[:,0], ene[:,[17, 18,19,20, 21, 22, 23, 24,25,26]])#, 9,10, 11, 12, 13, 14, 15, 16]])
plt.plot(fin[:,0], fin[:,1:])

#plt.savefig('plot_ene.png', dpi = 300)

print 'n. of points: ' + str(len(ene))
plt.show()
# In[15]:


# Evaluate expected values

plt.plot(expected_df.index, expected_df.iloc[:,1:])
#plt.plot(expected_df.index, expected_df.iloc[:,[8,9,10,]])
#plt.plot(expected_df.index, expected_df.iloc[:,[1,2,3,4,5,6,7]])
#plt.savefig('plot_prev.png', dpi = 300)

#print 'n. of points: ' + str(len(ene_prev))
plt.show()
#print expected_df


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


# In[3]:


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
    #expec_out.write(str(v) + '\t')



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



def extract_poly(i1, indexes, poly_cont):
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
    k1 = float(k1)
    k2 = float(k2)
    H = float(H)
    x = k1 + 2*k2*H
    return x
# In[7]:

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
    
    poli_out.write(lvls_list[i0] + '  k2 k1 k0\n' )
    
    
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
        pos = order_df.iloc[i1, i0]
        if pos == 'nan':
            five_E.append(float(expected_df.iloc[i1, i0]))
        
        else:
            pos = int(pos)
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
            vec = ['{:.7f}'.format(i) for i in vec_poly]
            poli_out.write('H: ' + str(H_value) + ' ' + vec[0] + ' ' + vec[1] + ' ' + vec[2] +'\n')
            #poli_out.write('H: ' + str(H_value) + ' ' + str(vec_poly) +'\n')
            
            
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
                order_df.iloc[i1+1, i0] = 'nan'#int(num[0])+1
                
                pass
            
            five_E.pop(0)
            five_H.pop(0)
    
    poli_out.write('\n')


# Write outputs
order_df.to_csv(path + 'order.txt', header = lvls_list, sep='\t', na_rep='na')
expected_df.to_csv(path + 'expected.txt', header = lvls_list, sep='\t', na_rep='na')
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


#Create the final_df array using the labels from order_df as a mask over ene_df

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
final_df.to_csv(path + 'final.txt', header = lvls_list, sep='\t', na_rep='na', float_format='%.7f' )


# In[]:
#
##
### Read poli.txt file
##
#

# Opens and reads
in_f = open(path + 'poli.out', 'r')
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

    
# In[]:
#
## Define variables
#
thrs_AE = float(0.5)
#thrs = float((raw_input('Valor max AE (thrs_AE): '))

thrs_pend = float(0.1)
#thrs_k1 = float((raw_input('Valor max d(pend) (thrs_pend): '))

out_f = open('transitions_sum.txt', 'w')

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
        sEs = final_df.iloc[i+4]
        H_value = H_values[i+4]
        H_middle = H_values[i+2]
        
        for a,b in combinations(sEs,2):  ### campo  middle
            AE_abs = abs(a-b)
            AE = a-b
            
            if AE_abs <= thrs_AE:
                
                #prints 'niveles buenos'
                index_a = sEs[sEs == a].index[0]
                index_b = sEs[sEs == b].index[0]
            
                # Saves levels somewhere, to then evaluate change in signs
                k1_a = float(dic[index_a][H_value][1])
                k1_b = float(dic[index_b][H_value][1])
                k2_a = float(dic[index_a][H_value][2])
                k2_b = float(dic[index_b][H_value][2])
                
                # Calculates pend
                pend_a = first_d(k1_a, k2_a, H_middle)
                pend_b = first_d(k1_b, k2_b, H_middle)
                print 'pend a' + str(pend_a)
                print 'pend b' + str(pend_b)
                
                # Pend differences
                Apend = abs(pend_a - pend_b)
                div_pend = pend_a / pend_b
                
                Abs_pend = abs(pend_a) - abs(pend_b) 
                
                ####### Signos diferentes dividir pends; if div < 0 entra
                               
                if pend_a > 0:
                    if pend_b < 0:
                        print 'estos molan'
                        if Apend <= thrs_pend:
                            # Stores
                            dic_2[H_value][index_a].append(index_b)
                            
                            # Writes
                            H_middle = '{:.7f}'.format(H_middle)
                            a = '{:.7f}'.format(a)
                            b = '{:.7f}'.format(b)
                            pend_a = '{:.7f}'.format(pend_a)
                            pend_b = '{:.7f}'.format(pend_b)
                            AE = '{:.7f}'.format(AE)
                            Apend = '{:.7f}'.format(Apend)
                    
                            out_f.write(H_middle + '     ' +  str(index_a) +'     ' +
                                        str(index_b)+ '     ' + a + '     ' + b+ '     ' +
                                        pend_a + '     ' + pend_b + '     ' + '     ' + AE+
                                        '     ' + Apend + '\n') 
                        
                elif pend_a < 0:
                    if pend_b > 0:
                        print 'estos molan'
                        if Apend <= thrs_pend:
                            # Stores
                            dic_2[H_value][index_a].append(index_b)
                            
                            # Writes
                            H_middle = '{:.7f}'.format(H_middle)
                            a = '{:.7f}'.format(a)
                            b = '{:.7f}'.format(b)
                            pend_a = '{:.7f}'.format(pend_a)
                            pend_b = '{:.7f}'.format(pend_b)
                            AE = '{:.7f}'.format(AE)
                            Apend = '{:.7f}'.format(Apend)
                            
                            out_f.write(H_middle + '     ' +  str(index_a) +'     ' +
                                        str(index_b)+ '     ' + a + '     ' + b+ '     ' +
                                        pend_a + '     ' + pend_b + '     ' + '     ' + AE+
                                        '     ' + Apend + '\n')
                

out_f.close()


# In[ ]:
#################
#    OLD
#################
                    
#
##
### Functions
##
#
def pend(l_real, l_imag, l_projections, n):
    g = float(7)/6
    beta = 0.46686
    l_mj = []
    
    for i in l_projections:
        vec = i.split('/')
        l_mj.append(float(vec[0]))
    l = list(range(0,n))
    coef_l = []
    
    for i in l:
        x = complex(float(l_real[i]), float(l_imag[i]))    # complex
        coef_i = g*beta*(abs(x)**2)*l_mj[i]  # complex
        coef_l.append(coef_i)
    coef = sum(coef_l)
    
    return coef


def compare_pend(d_pends, d_eigenvectors): #d_eigenvalues
    d_dif = {}
    d_uniq = {}
    lim = 10**(-3)
    
        
    for a, b in combinations(d_pends, 2):
        #difE != 0
        # continue
        if d_eigenvectors[a] != d_eigenvectors[b]:
            value = abs(d_pends[a]- d_pends[b])

            if value < lim:
                d_dif[a + b] = value
    
    return d_dif


def sderivative(pend1, pend2, field1, field2):
    value = (float(pend1) - float(pend2)) / (float(campo1) - float(pend2))



# In[98]:


#
##
### Slope calculation module from eigenvalues and eigenvectors
##
#

# Out files
out_f1 = open(path + 'global.out', 'w')
out_f3 = open(path + 'pends.out', 'w')
out_f3.write('H_value\t')

for i in list(range(1,tot_lvls+1)):
    out_f3.write('lvl_' + str(i) + '\t')
out_f3.write('\n')


#
## Parsing input file 'simpre.out'
#
with open(path + 'simpre.out', 'r') as f:
    lines = f.read().split('\n')




d_H = {}
H_values = []
d_eivec_pos = {}


c = 0
for i in lines:
    c+=1
    if 'Magnetic Field' in i:
        vec = i.split()
        H = vec[4]
        H_values.append(H)
        d_H[H] = c
        
    elif 'Eigenvectors' in i:
        d_eivec_pos[H] = c


c2 = 0
for i in H_values:
    c2 += 1
    out_f1.write('#################################################\n')
    out_f1.write('## Magnetic Field value = ' + str(i) + '\n')
    out_f1.write('#################################################\n\n')
    out_f1.write('## tag   EigValues   SpinProjections   pend   EigVect(R I)\n')
    
    out_f2 = open(path + 'temp_' + str(c2) + '.out', 'w')
    out_f2.write('#################################################\n')
    out_f2.write('## Magnetic Field value = ' + str(i) + '\n')
    out_f2.write('#################################################\n\n')
    out_f2.write('## tag   EigValues   SpinProjections   pend   EigVect(R I)\n')
    
    out_f3.write(str(i) + '\t')
    
    
    d_eigenvalues = {}
    d_eigenvectors = {}
    d_pend = {}
    l_pend = []
    
    eigenvalues_list = lines[(d_H[i] + 18) : (d_H[i] + 44) ] # 18 y 44 son valores fijos 
    
    for i2 in eigenvalues_list:
        vec = i2.split()
        d_eigenvalues['lvl_' + str(vec[3])] = vec[0]
        
    eigenvectors_list = lines[d_eivec_pos[i] : d_eivec_pos[i] + tot_lvls] # tot_lvls valor variable
    
    c3 = 0
    for i3 in eigenvectors_list:
        c3 += 1
        d_eigenvectors['lvl_' + str(c3)] = i3
        
    
    
    #
    ## Calculates slopes
    #
    for i4 in d_eigenvectors:
        vec = d_eigenvectors[i4].split()
        c_real =  vec[0::2]
        c_imag = vec[1::2]

        pend_i = pend(c_real, c_imag, projections_l,tot_lvls)
        d_pend[str(i4)] = pend_i
        l_pend.append(pend_i)
        

    ## Condition 1 : discriminates degeneracy and calculates difference between pends
    d_diffs = compare_pend(d_pend, d_eigenvalues)
    
    
    
    #### Writing temp_f and global .out
    for i in list(range(1,tot_lvls+1)):
        tag = 'lvl_' + str(i)
        out_f1.write(tag + ' ' + str(d_eigenvalues[tag]) + ' ' +  str(projections_l[i-1])
                    + ' ' + str(d_pend[tag]) + ' ' + str(d_eigenvectors[tag]) +  '\n')  # Falta Col [0/1] compare pend
        out_f2.write(tag + ' ' + str(d_eigenvalues[tag]) + ' ' +  str(projections_l[i-1])
                    + ' ' + str(d_pend[tag]) + ' ' + str(d_eigenvectors[tag]) +  '\n')  # Falta Col [0/1] compare pend
        
        out_f3.write(str(d_pend[tag]) + '\t')
    
    
    
    out_f1.write('\n\n\n')
    out_f1.write('## Differences map\n')
    out_f1.write('\n\n\n')
    
    out_f2.write('\n\n\n')
    out_f2.write('## Differences map\n')
    
    for i in d_diffs:
        out_f1.write(i + ' ' + str(d_diffs[i]) + '\n')
        out_f2.write(i + ' ' + str(d_diffs[i]) + '\n')
    
    out_f2.close()
    
    out_f3.write('\n')
    
    
    
    ## Break
    if c2 == 47:
        break



out_f1.close()
out_f3.close()


