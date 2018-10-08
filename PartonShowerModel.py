# -*- coding: utf-8 -*-
"""
Created on Mon May  7 21:52:54 2018

@author: outong
"""
import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt

partontype=np.dtype({'names':['energy','xf','yf','ang'],\
                     'formats':['f','f','f','f']})

norm_z=1/np.log(2)
'''
Generate z distributed as PDF(z)=1/(z+1)
 with Inverse Transform Method
'''
def gen_z():
    y=random.random()
    z=np.exp(y/norm_z)-1
    return z

norm_theta=1/np.log(1+np.pi/2)
'''
Generate theta distributed as PDF(theta)=1/(theta+1)
 with Inverse Transform Method
'''
def gen_theta():
    y=random.random()
    t=np.exp(y/norm_theta)-1
    return t

'''
Drift time of the initial parton
'''
def dt(z,ang):
    return (z+1)*(ang+1)


E0=100.0 #(Gev) Energy of the initial parton (just an assumption)
Ecri=20.0 #(Gev) Energy threshold
x0=0.0 #The position of the initial parton
y0=0.0
theta0=0.0 
A=np.array([(E0,x0,y0,theta0)],dtype=partontype)
print(A.dtype)
i=0

while 1:
    E0=A[i]['energy']
    
    #If E<Ecri, stop the splitting process
    if (i==len(A)-1)and(E0<Ecri):
        X=[A[i]['xf']]
        Y=[A[i]['yf']]
        
        xf=A[i]['xf']+1
        yf=A[i]['yf']+1*A[i]['ang'] 
        
        X.append(xf)
        Y.append(yf)
        
        plt.plot(X,Y)
       # plt.scatter(X,Y)
        del X,Y
        break
    elif (E0<Ecri):
        X=[A[i]['xf']]
        Y=[A[i]['yf']]
        
        xf=A[i]['xf']+1
        yf=A[i]['yf']+1*A[i]['ang'] 
        
        X.append(xf)
        Y.append(yf)
        
        plt.plot(X,Y)
        del X,Y
    
        i+=1
        continue
    
    ang0=A[i]['ang']
    
    #generate a pair of [z,theta]
    z=gen_z()
    theta=gen_theta()
    
    beta=1 #Relativistic limit
    xf=beta*dt(z,ang0)#small angle approximation
    yf=beta*dt(z,ang0)*ang0
    
    X=[A[i]['xf']]
    Y=[A[i]['yf']]
    
    A[i]['xf']+=xf
    A[i]['yf']+=yf
    
    X.append(A[i]['xf'])
    Y.append(A[i]['yf'])
    
    plt.plot(X,Y)
    #plt.scatter(X,Y)
    del X,Y
    
    #Determine the track of the generated partons 
    #according to conservation of energy and momenta
    E1=(1-z*1/2)*E0 #considering E1 is the more energetic one
    E2=E0-E1
    ang1=ang0+theta
    ang2=ang0-E1*theta/E2
    
    B=np.array([(E1,A[i]['xf'],A[i]['yf'],ang1)],dtype=partontype)
    C=np.array([(E2,A[i]['xf'],A[i]['yf'],ang2)],dtype=partontype)
    A=np.concatenate((A,B,C))

    i+=1

plt.xlabel("X",fontsize=15)
plt.ylabel("Y",fontsize=15)
plt.show()
       