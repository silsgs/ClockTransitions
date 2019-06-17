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
out1_file = path + (sys.argv[1].replace('.txt','.fderivative'))
out2_file = path + (sys.argv[1].replace('.txt','.sderivative'))

total_lines = os.system('cat ' + in_file + ' | wc -l')

#---------------------------------#
## Functions ##

# plot energy levels
#def ene_plot(dct):


# Calculates the first derivative
def derivative(n, list1, list2):
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

d_dct = {}
for i in dct:
    d = derivative(total_lines, field_vls, dct[i])
    d_dct['d_%s' % i] = d

sd_dct = {}
for i in d_dct:
    sd = derivative(total_lines-1, field_vls, d_dct[i])
    sd_dct['sd_%s' % i] = sd

#print sd_dct

# Writes output1 - first derivative
out1_f = open(out1_file, 'w')
out1_f.write('%f   %s\n' % (float(field_vls[0]), str(' - ')))


for n in list(range(len(d_dct['d_lv_1']))):
    row1 = [field_vls[n+1] + '   ']
    for i in list(range(1, len(d_dct)+1)):
        row1.append(str(d_dct['d_lv_%s' % i][int(n)]))
    for item in row1:
        out1_f.write('%f  ' % float(item))
    out1_f.write('\n')



# Writes output2 - second derivative
out2_f = open(out2_file, 'w')
out2_f.write('%f   %s\n' % (float(field_vls[0]), str(' - ')))
out2_f.write('%f   %s\n' % (float(field_vls[1]), str(' - ')))

print len(sd_dct)

c = 0
for n in list(range(len(sd_dct['sd_d_lv_1']))):
    row2 = [field_vls[n+2] + '   ']
    for i in list(range(1, len(sd_dct)+1)):
        row2.append(str(sd_dct['sd_d_lv_%s' % i][int(n)]))
    for item in row2:
        out2_f.write('%f  ' % float(item))
    out2_f.write('\n')


out1_f.close()
out2_f.close()
