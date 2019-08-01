import sys
import os
#import numpy as np
#import pyplot

#---------------------------------#
# Parameters defined by the user
#levels = [1,2]
total_levels = 5
#total_lines = 191

#---------------------------------#
## Definitions
path = os.getcwd() + '/'
in_file = path + sys.argv[1]
out_file = path + (sys.argv[1].replace('.txt','.derivatives'))

total_lines = os.system('cat ' + in_file + ' | wc -l')

#---------------------------------#
## Functions ##

# plot energy levels
#def ene_plot(dct):


# Calculates the first derivative
def first_d(n, list1, list2):
    d = []
    for i in list(range(0, n)):
        if i == 0:
            continue
        else:
            di = (float(list2[i]) - float(list2[i-1]))/(float(list1[i]) - float(list1[i-1]))
            d.append(di)
    return d



#---------------------------------#
#Read and parse input file
in_f = open(in_file, 'r')

total_lines = 0
field_vls = []
dct = {}
for i in list(range(1, total_levels)):
    dct['lv_%s' % i] = []
#print dct

for line in in_f:
    vec = line.strip().split()
    #total_levels = len(vec) - 1
    field_v = vec[0]
    field_vls.append(field_v)
    total_lines += 1

    for i in list(range(1, total_levels)):
        dct['lv_%s' % i].append(vec[i])

in_f.close()




#---------------------------------#
## Sets of data
#dct{'lv_1':[0,..,n],'lv_2':[0,..,n],....,'lv_n':[0,..,n]}
#field_vls=[]
#d_dct{}




#---------------------------------#
#Function calls and out writing

d_dct = {}
for i in dct:
    d = first_d(total_lines, field_vls, dct[i])
    d_dct['d_%s' % i] = d

#print d_dct
#print dct['lv_1'][0]

out_f = open(out_file, 'w')
out_f.write('%f   %s\n' % (float(field_vls[0]), str(' - ')))


for n in list(range(len(d_dct['d_lv_1']))):
    row = [field_vls[n+1] + '   ']
    for i in d_dct:
        row.append(str(d_dct[i][int(n)]))
    for item in row:
        out_f.write('%f  ' % float(item))
    out_f.write('\n')

out_f.close()
