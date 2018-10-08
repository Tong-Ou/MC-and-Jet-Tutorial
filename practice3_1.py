# -*- coding: utf-8 -*-
"""
Created on Tue May  1 11:04:35 2018

@author: Tong Ou
"""

import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Gauss function (for fitting)
def gauss(x,norm,mu,sigma):
    t=norm*np.exp(-(x-mu)**2/(2*sigma**2))
    return t

n_points=1000
nr=1000
a=np.zeros(1000)
for i in range(nr):
    nc=0
    for j in range(n_points):
        x=random.random()-1/2
        y=random.random()-1/2
        #generate a point in the square of side=1
        if x**2+y**2<=1/4:
            nc+=1
            #pick out the points in the circle of diameter=1
    pi=4*nc/n_points
    a[i]=pi

'''
Repeat the calculation of pi for nr times,
draw all the results to a histogram and fit it with a Gauss
function.
'''
n_bins=50
n,bins,patch=plt.hist(a,bins=n_bins,facecolor="green")
popt,pcov=curve_fit(gauss,bins[1:],n)
mu=popt[1]
sigma=abs(popt[2])
mu2=str(round(mu,3))
sigma2=str(round(sigma,4))
plt.plot(bins[1:],gauss(bins[1:],*popt),'r',\
         label='$\mu=$'+mu2+'\t'+'$\sigma=$'+sigma2)
plt.xlabel("Calculated Pi",fontsize=15)
plt.ylabel("Counts",fontsize=15)
plt.legend()
plt.show()
print(mu,sigma)
