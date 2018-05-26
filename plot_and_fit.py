# -*- coding: utf-8 -*-
"""
Created on Tue May 22 21:05:34 2018

@author: outong
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def gauss(x,norm,mu,sigma):
    t=norm*np.exp(-(x-mu)**2/(2*sigma**2))
    return t

a=np.loadtxt('njet_p-1_R1.0.txt')
nn=len(a)
nj=[]
pm=[]
for i in range(nn):
    nj.append(a[i,0])
    pm.append(a[i,1])
    
fig,(ax0,ax1)=plt.subplots(ncols=2,figsize=(14,5))
    
n_bins=50
n,bins,patch=ax0.hist(nj,bins=n_bins,facecolor="purple")


p0=[150,3,1]
popt,pcov=curve_fit(gauss,bins[1:],n,p0)
norm=popt[0]
mu=popt[1]
sigma=abs(popt[2])
mu2=str(round(mu,2))
sigma2=str(round(sigma,3))
ax0.plot(bins[1:],gauss(bins[1:],*popt),'r',\
         label='$\mu=$'+mu2+'\t'+'$\sigma=$'+sigma2)

ax0.set_xlabel("Number of Jets",fontsize=15)
ax0.set_ylabel("Counts",fontsize=15)
ax0.legend()


n2,bins2,patch2=ax1.hist(pm,bins=n_bins,fill=False,color='g')
p2=[100,1,1]
#popt2,pcov2=curve_fit(gauss,bins2[1:],n2,p2)
#mu_pm=str(round(popt[1],2))
#sigma_pm=str(round(popt[2],4))
#ax1.plot(bins2[1:],gauss(bins2[1:],*popt2),'r',\
#         label='$\mu=$'+mu_pm+'\t'+'$\sigma=$'+sigma_pm)
ax1.set_xlabel("Pseudomass of Jets",fontsize=15)
ax1.set_ylabel("Counts",fontsize=15)

ax1.legend()

plt.show()
print(norm,mu,sigma)