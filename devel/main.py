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


# =============================================================================
#                              FLOW
# =============================================================================
# 1. Run your simpre calculation and plot your results.
# 2. Define a trainning data set region without non-avoided crossings
#     - Minimum of 3-5 points (plot your results)
# 3. Run Poly.py 
# 4. Run expected_values.py
# 5. Run main.py


# In[15]:
#
## Functions module
#
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
    #expec_out.write(str(v) + '\t')


def truncate(f, n):
    '''Truncates/pads a float f to n decimal places without rounding'''
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return '.'.join([i, (d+'0'*n)[:n]])


def search(v0, v1, v2, list_e): #index = pos-1
    '''Finds the corresponding E value in simpre.ene comparing with expected_value'''
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

# In[4]:
# 1. Plot your results
ene_f = path + "simpre.ene"
ene = np.loadtxt(ene_f, dtype= float)

# Check if subdir_path exists and if its empty
plotsdir_path = path + 'plots/'
if not os.path.exists(plotsdir_path):
    os.makedirs(plotsdir_path)

plt.plot(ene[:,0], ene[:,1:])
print 'n. of points: ' + str(len(ene))
plt.xlabel(u'Magnetic field (T)')
plt.ylabel(u'Energy (cm-1)')
plt.savefig('plots/plot_ene.png', dpi = 300)
plt.show()

#
## Initial definitions
#
# Defined by the user so far; automatize for the future
#J = float(8)
#I = float(7.5)
J = float(raw_input('Valor de J: '))
I = float(raw_input('Valor de I: '))
g_par = 0.0000
g_per = 0.0000



## Lists
# Total spin levels of the system
tot_lvls = int((2*J + 1)*(2*I + 1))
m_j = list(range(int(-(2*J)), int((2*J+2)), int(2)))
m_i = list(range(int(-(2*I)), int((2*I+2)), int(2)))



# Loads new 'simpre.ene' file with energies of each level (quantum number) at diff H field
n_levels = raw_input('How many n lower-energy levels are you interested in?:')

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



## Generate projections list from m_j and m_i
projections_l = []
for j in m_j:
    for i in m_i:
        proj = str(j) + '/' + str(i)
        projections_l.append(proj)



   
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
        '''
        print 'valor de campo actual'
        print H_value
        print 'valor de campo next'
        print H_next
        '''
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
            '''
            print 'aqui los 5 campos'
            print five_H
            print 'aqui las 5 Es'
            print five_E
            '''
            
            # Fitting
            vec_poly = poly(five_H, five_E)
            p = np.poly1d(vec_poly)
            d_poli[i1] = p
            vec = ['{:.7f}'.format(i) for i in vec_poly]
            poli_out.write('H: ' + str(H_value) + ' ' + vec[0] + ' ' + vec[1] + ' ' + vec[2] +'\n')
            
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
final_df.to_csv(path + 'final.txt', header = lvls_list, sep='\t', na_rep='na', float_format='%.7f' )

# In[]:
# Visual checking of results final vs. expected
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey = True)
ax1.plot(ene[:,0], ene_df.iloc[:,0:])
ax1.set_title('Energies')
ax2.plot(expected_df.index, expected_df.iloc[:,1:])
ax2.set_title('Expected energies')
ax3.plot(final_df.index, final_df.iloc[:,1:])
ax3.set_title('Final energies')
fig.text(0.5,0.01, 'Magnetic field (T)', ha='center', fontsize = 12)

plt.savefig('plots/summary_plots.png', dpi = 300)

# In[]:
#
##
### Read amd parse poli.txt file
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

# Open output files
out_avoid = open('avoided_sum.txt', 'w')
out_nonavoid = open('nonavoided_sum.txt', 'w')


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
        
        for a,b in combinations(sEs,2):  ### campo  middle
            # Calculates AE
            AE = abs(a-b)
            
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
res = open(path + 'results.out', 'w')

# Write header res file
res.write('#Level_a     Level_b     Type_crossing     H_crossing     E_crossing1     E_crossing2     Curvature\n')

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
        
        # 
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
                H_cross = '{:.5f}'.format(H_cross)
                E_cross1 = '{:.5f}'.format(E_cross1)
                E_cross2 = '{:.5f}'.format(E_cross2)
                curvature = calc_curvature(H_cross, k1_a, k1_b, k2_a, k2_b)
                curv1 = '{:.5f}'.format(curvature[0])
                curv2 = '{:.5f}'.format(curvature[1])
                #curvature = '{:.5f}'.format(curvature)
                res.write(a + '     '+ b + '     '+ type_cross + '     '+ H_cross + '     '+
                          E_cross1 + '     ' + E_cross2 + '     ' + curv1 +'     ' + curv2 +  '\n' )
                
                
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
                H_cross = '{:.5f}'.format(H_cross)
                E_cross1 = '{:.5f}'.format(E_cross1)
                E_cross2 = '{:.5f}'.format(E_cross2)
                curvature = calc_curvature(H_cross, k1_a, k1_b, k2_a, k2_b)
                curv1 = '{:.5f}'.format(curvature[0])
                curv2 = '{:.5f}'.format(curvature[1])
                res.write(a + '     '+ b + '     '+ type_cross + '     '+ H_cross + '     '+
                          E_cross1 + '     ' + E_cross2 + '     ' + curv1 +'     ' + curv2 +  '\n' )
                

res.close()