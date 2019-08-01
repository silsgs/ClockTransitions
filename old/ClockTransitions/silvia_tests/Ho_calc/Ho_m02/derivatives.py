import sys
import os
#import numpy as np
#import pyplot

#---------------------------------#
# Parameters defined by the user
levels = [6,12]
total_levels = 17
#total_lines = 191

#---------------------------------#
## Definitions
path = os.getcwd() + '/'
in_file = path + sys.argv[1]
out_file = path + (sys.argv[1].replace('.ene','.derivatives'))

total_lines = os.system('cat ' + in_file + ' | wc -l')

#---------------------------------#
## Functions ##


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
#Function calls and out writing

out_f = open(out_file, 'w')

for i in levels:
    if len(levels) < 1:
        print 'Error: you should indicate the levels'

    elif len(levels) == 1:
        d = first_d(total_lines, field_vls, dct['lv_%s' % i])
        dct['d_%s' % i] = d
        out_f.write(str(dct['d_%s' % i]))

    else:
        d = first_d(total_lines, field_vls, dct['lv_%s' % i])
        dct['d_%s' % i] = d
        print d

    for n in list(range(len(d))):
        out_f.write('%f   %f\n' % (float(field_vls[n+1]), float(d[n])))

out_f.close()
