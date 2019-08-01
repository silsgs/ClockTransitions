import sys
import os
import numpy as np
import matplotlib.pyplot as plt

path = os.getcwd() + '/'
in_file = path + 'simpre.ene'
#print in_file
in_f = open(in_file, 'r')

#plt.figure(1)

data = np.loadtxt(in_file)    
#plt.subplot(311)
plt.plot(data)
plt.show(block=True)

#plt.figure(2)
derivatives_f = path + 'simpre_ene_test1.fderivative'
data2 = np.loadtxt(derivatives_f, dtype= float, skiprows = 1)    
#plt.subplot(312)
plt.plot(data2)
plt.show(block=True)

#plt.figure(3)
derivatives_f = path + 'simpre_ene_test1.sderivative'
data3 = np.loadtxt(derivatives_f, dtype= float, skiprows = 2, usecols = (1,2,3,4))    
#plt.subplot(313)
plt.plot(data3)
plt.show(block=True)
