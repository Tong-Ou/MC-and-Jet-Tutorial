# -*- coding: utf-8 -*-
"""
Created on Wed May  2 19:51:33 2018

@author: outong

Generate random numbers distributed as exp(-x)
with Inverse Transform Method.
"""
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt

#define the CDF of distribution exp(-x)
def cdf(x,xmin):
    t=N*(np.exp(-xmin)-np.exp(-x))
    return t

#define the inverse transform of CDF of exp(-x)
def re_cdf(y,xmin):
    t=-np.log(np.exp(-xmin)-y/N)
    return t

xmin=0
xmax=5
N=1/(np.exp(-xmin)-np.exp(-xmax))#normalization constant
n_points=1000
a=np.zeros(n_points)
for i in range(n_points):
    y=random.random()
    x=re_cdf(y,xmin)#x is the inverse transform of y
    a[i]=x
n_bins=50
plt.hist(a,bins=n_bins,normed=1)
x=np.linspace(0,5,100)
pdf=N*np.exp(-x)

#drawing the PDF curve for reference.
plt.plot(x,pdf,'r',label="PDF(x)=$N\cdot e^{-x}$")
plt.xlabel("X",fontsize=15)
plt.ylabel("PDF(X)",fontsize=15)
plt.legend()
plt.show()
    
