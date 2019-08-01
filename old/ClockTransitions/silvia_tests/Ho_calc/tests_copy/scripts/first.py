### Extracting R_values and I_values 
## from simpre.out for XXX levels w/ lower E
#
J = float(6)
I = float(0.5)
#
##

tot = int((2*J + 1)*(2*I + 1))

c = 0
list_eigen = []

out_f = open(path + 'short.out', 'w')
with open(path + 'simpre.out', 'r') as f:
    lines = f.read().splitlines()
    last_lines = lines[-11:]


## Create d_eigenvalues    
for i in lines:
    if 'Eigenvalues' in i:
        c+=1
        #print 'True1'
        if c > 1 :
            #print 'break1'
            break
        if line.startswith('Basis Functions'):
            #print 'break2'
            break
        else:
            ind = lines.index(i)
            #print 'else1  ', ind

list_eigen = lines[(ind + 2):(tot + ind + 2)]
d_eigen = {}

for i in list_eigen:
    vec = i.split()
    d_eigen[vec[3]] = vec[0]


# Create d_real & d_imag
c=0
d_real = {}
d_imag = {}
for i in last_lines:
    vec = i.split()
    c += 1
    d_real[c] = vec[0::2]
    d_imag[c]= vec[1::2]


lc = list(range(1,c+1))
m_j = list(range(-12,14,2))
m_i = list(range(-1,2,2))


## Output writing
out_f.write('## total levels:' + str(len(d_real)) + '\n')

for l in lc:
    for e in d_real[l]:
        out_f.write(str(e) + ' ')
#    out_f.write('\n')
    
    for e_i in d_imag[l]:
        out_f.write(str(e_i) + ' ')
    out_f.write('\n')
    
    for j in m_j:
        mj = float(j)/2
        if len(m_i)>0:
            for i in m_i:
                mi = float(i)/2
                out_f.write(str(mj) + '/' + str(mi) + '  ')
        else:
            out_f.write(str(mj) + '/0  ')
    out_f.write('\n\n')

out_f.close()
