## Conditions
#

# First condition : levels with unique eigenvalue 
d_uniq = {}

for i in list(range(1,len(d_eigen))):
    d_eigen['0'] = 0.0
    value = float(d_eigen[str(i)]) - float(d_eigen[str(i-1)])
    
    if value != 0.0: # if AE < 18*-8: mismo nivel 
        d_uniq['lvl_' + str(i)] =  d_eigen[str(i)]
    
    if float(d_eigen['1']) == float(0.0): 
        d_uniq['lvl_1'] = d_eigen['1']


print d_uniq

print d_pend
#print d_eigen['1']
#print float(0.0)


# Second condition: 
# - recorrer valores de campo de los niveles unicos
# - calcular pendiente en cada uno de ellos (funcion 'coef')
# - calcular diferencia entre ppendientes en cada punto


# pend_i = pend_j if Apend < AE/AH; 10**-8-10**-4





###########################################################################
10/04/2019
###########################################################################



##### Run at the same folder as simpre calculation

## Defined by the user so far; automatize for the future
J = float(6)
I = float(0.5)
tot_lvls = int((2*J + 1)*(2*I + 1))
m_j = list(range(int(-(2*J)), int((2*J+2)), int(2)))
m_i = list(range(int(-(2*I)), int((2*I+2)), int(2)))

## Parse .par
with open(path + 'simpre.par', 'r') as par_f:
    #lines = f.read().split('\n')
    for i in par_f:
        if 'fieldstart' in i:
            fieldstart = i[19 + len('fieldstart'):18 + len('fieldstart')+ 7]
            #print fieldstart
        elif 'fieldend' in i:
            fieldend = i[19 + len('fieldend'):18 + len('fieldend')+ 7]
            #print fieldend
        elif 'fieldstep' in i:
            fieldstep = i[19 + len('fieldstep'):18 + len('fieldstep')+ 8]
            #print fieldstep
par_f.close()


## Parse 'simpre.out'
# Objects: para cada valor del barrido campo magnetico (de 'simpre.out')
# d_eigenvalues ; d_eigenvectors ; d_real ; d_imag ; d_pend ; d_uniq


#
## Parsing input file 'simpre.out'
#

out_f = open(path + 'temp1.out', 'w')
with open(path + 'simpre.out', 'r') as f:
    lines = f.read().split('\n')

c = 0
d_H = {}
H_values = []
d_eivec_pos = {}

for i in lines:
    c+=1
    if 'Magnetic Field' in i:
        vec = i.split()
        H = vec[4]
        H_values.append(H)
        d_H[H] = c
        out_f.write(i + '\n')
        
    elif 'Eigenvectors' in i:
        d_eivec_pos[H] = c
#print d_eivec_pos    


### Main loop
#n = list(range(len(H_values)))
#print H_values
c2 = 0

for i in H_values:
    d_eigenvalues = {}
    d_eigenvectors = {}
    
    c2 += 1
    
    eigenvalues_list = lines[(d_H[i]+18) : (d_H[i] + 44) ] # 18 y 44 son valores fijos 
    
    for i2 in eigenvalues_list:
        vec = i2.split()
        d_eigenvalues['lvl_' + str(vec[3])] = vec[0]
        
    #print lines[d_eivec_pos[i] : d_eivec_pos[i] + 26]
    
    eigenvectors_list = lines[d_eivec_pos[i] : d_eivec_pos[i] + tot_lvls] # tot_lvls valor variable
    #print eigenvectors_list
    
    c3 = 0
    for i3 in eigenvectors_list:
        c3 += 1
        d_eigenvectors['lvl_' + str(c3)] = i3
        d_pend = {}
        
        for i4 in d_eigenvectors:
            #print len(d_eigenvectors[i4].split())
            vec = d_eigenvectors[i4].split()
            c_real =  vec[0::2]
            c_imag = vec[1::2]
            #print 'c_real:', c_real
            #print 'c_imag:', c_imag
            
            d_pend[i3] = pend(c_real, c_imag, tot_lvls)
            print d_pend
            
        
    
    
    
    #print d_eigenvectors['lvl_1']
    
    if c2 == 1:
        break
"""
    for i3 in eigenvectors_list:
        vec



    


#print d_eigenvalues        
#print eigenvalues_list
#print lines[(d_H['-0.34500'] + 18):(d_H['-0.34500'] + 44)]
#for i in 
"""
out_f.close()

