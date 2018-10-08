# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 20:02:12 2018

@author: Tong Ou
"""
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt

n_points=1000

def gauss(x,mu,sigma):
    t=1/np.sqrt(2*np.pi*sigma**2)*np.exp(-(x-mu)**2/(2*sigma**2))
    return t

#Generate random numbers of Gauss distribution
def gen(xmin,xmax,sigma,mu):
    a=[]
    L=xmax-xmin
    i=0
    while i in range(n_points):
        t=random.random()
        x=xmin+L*t
        y=random.random()
        s=gauss(x,mu,sigma)/gauss(mu,mu,sigma)
        if y<=s:
            a.append(x)
            i+=1
    return a

a=gen(0,10,2,5)
n_bins=50
#drawing the random numbers to a histogram
plt.hist(a,bins=n_bins,normed=1)
plt.xlabel("X",fontsize=15)
plt.ylabel("PDF(X)",fontsize=15)
plt.title("bins=50",fontsize=15)

#drawing the Gauss function for reference
x=np.linspace(0,10,n_bins)
y=gauss(x,5,2)
plt.plot(x,y,'r',label="Gauss distribution")
plt.legend()
plt.show()

        

