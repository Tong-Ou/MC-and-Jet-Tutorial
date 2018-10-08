# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy.random as random
import matplotlib.pyplot as plt

n=1000 #give the number of required random numbers

#The function to generate n random numbers ranging in [0,1]
def gen_1():
    a=[]
    for i in range(0,n):
        x=random.random()
        a.append(x)
    return a

#The function to generate n random numbers ranging in [xmin,xmax]
def gen_2(xmin,xmax):
    b=[]
    L=xmax-xmin
    for i in range(0,n):
        t=random.random()
        x=xmin+L*t
        b.append(x)
    return b

#generate the random numbers
a=gen_1()
b=gen_2(5,15)

#drawing the number to histograms
n_bins=50
plt.hist(a,bins=n_bins)
plt.ylabel("Counts",fontsize=15)
plt.xlabel("X",fontsize=15)
plt.title("bins=50",fontsize=15)
plt.show()
plt.hist(b,bins=n_bins)
plt.ylabel("Counts",fontsize=15)
plt.xlabel("X",fontsize=15)
plt.title("bins=50",fontsize=15)
plt.show()

