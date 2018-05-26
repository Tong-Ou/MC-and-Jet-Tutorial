# -*- coding: utf-8 -*-
"""
Created on Wed May  2 21:42:30 2018

@author: outong

Generate random numbers distributed as exp(-x)
with accept-reject method.
"""

import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt

#accept-reject method
def ac_re(n_points):
    nc=0
    i=0
    ymax=N*np.exp(-xmin)
    a=np.zeros(n_points)
    while i in range(n_points):
        s=random.random()
        x=xmin+s*(xmax-xmin)
        y=ymax*random.random()
        if y<=N*np.exp(-x):
            a[i]=x
            i+=1
        nc+=1
    eff=n_points/nc
    return eff
        
'''
Investigate the relationship between 
accepted x values and efficiency.
'''
xmin=0
xmax=5
N=1/(np.exp(-xmin)-np.exp(-xmax))#normalization constant
NP=np.arange(50,1100,50)
L=len(NP)
EFF=np.zeros(L)
for i in range(L):
    EFF[i]=ac_re(NP[i])

plt.plot(NP,EFF)
plt.xlabel('Accepted X values',fontsize=15)
plt.ylabel('Efficiency',fontsize=15)
plt.show()

