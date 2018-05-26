# -*- coding: utf-8 -*-
"""
Created on Tue May  1 22:46:16 2018

@author: Tong Ou
"""

import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Gauss function (for fitting)
def gauss(x,mu,sigma):
    t=1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))
    return t

#calculating pi for nr times and get the mu and sigma
def calc_pi(n_points,nr):
    a=np.linspace(0,0,nr)
    for i in range(nr):
        nc=0
        for k in range(n_points):
            x=random.random()-1/2
            y=random.random()-1/2
            #generate a point in the square of side=1
            if x**2+y**2<=1/4:
                nc+=1
                #pick out the points in the circle of diameter=1
        pi=4*nc/n_points
        a[i]=pi
    n_bins=50
    n,bins,patch=plt.hist(a,bins=n_bins,normed=1,facecolor="green")
    popt,pcov=curve_fit(gauss,bins[1:],n)
    mu=popt[0]
    sigma=abs(popt[1])
    return mu,sigma

'''
Fix the number of random points=1000 and alter 
the calculation times to see how sigma and mu
evolve.
'''

n_points=1000
NR=np.arange(50,3100,200)
l=len(NR)
MU=np.zeros(l)
SIGMA=np.zeros(l)
for i in range(l):
    MU[i],SIGMA[i]=calc_pi(n_points,NR[i])
fig,ax1=plt.subplots()
ax1.set_xlabel("Calculation Times",fontsize=15)
ax1.set_ylabel("$\sigma$",fontsize=15)
ax1.plot(NR,SIGMA,label="$\sigma$")
ax1.set_ylim(0.0450,0.065)
ax1.set_xlim(0,3100)
ax1.legend()
ax2=ax1.twinx()
ax2.set_ylabel("$\mu$",fontsize=15)
ax2.plot(NR,MU,'r',label="$\mu$")
ax2.legend()

