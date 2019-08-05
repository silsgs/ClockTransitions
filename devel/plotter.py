#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 14:24:18 2019

@author: silvia
"""


import os
import numpy as np
import matplotlib.pyplot as plt


path = os.getcwd() + '/'

ene_f = path + "simpre.ene"
ene = np.loadtxt(ene_f, dtype= float)



if __name__ == "__main__":
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
