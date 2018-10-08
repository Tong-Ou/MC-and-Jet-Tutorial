# -*- coding: utf-8 -*-
"""
Created on Sat May 19 10:05:39 2018

@author: outong
"""

import numpy as np
import numpy.random as random
import matplotlib.pyplot as plt

partontype=np.dtype({'names':['energy','px','py','pz','theta','phi'],\
                     'formats':['f','f','f','f','f','f']})

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
Coordinate transform (from natural coordinates to lab frame)
'''
def R(X,theta,phi):
    c1=np.cos(theta)
    c2=np.cos(phi)
    s1=np.sin(theta)
    s2=np.sin(phi)
    RY=np.array([[1,0,0,0],[0,c1,0,s1],[0,0,1,0],[0,-s1,0,c1]])
    RZ=np.array([[1,0,0,0],[0,c2,-s2,0],[0,s2,c2,0],[0,0,0,1]])
    RT=np.dot(RY,RZ)
    Xnew=np.dot(RT,X)
    return Xnew
        
'''
Generate the energy of one initial parton according to
the distribution of exp(-alpha*E)
'''
def gen_E0():
    alpha=0.001 #an assumption
    Emax=1000 #(GeV)
    norm=alpha/(1-np.exp(-alpha*Emax))
    y=random.random()
    E=-1/alpha*np.log(1-alpha*y/norm)
    return E

'''
3D Parton Shower generation
'''
def PartonShower3D(A,file,s):
    f=open(file,s)
    i=0
    while 1:
        E0=A[i]['energy']
        if (i==len(A)-1)and(E0<Ecri):
            print(E0)
            E0=round(E0,4)
            px=round(A[i]['px'],4)
            py=round(A[i]['py'],4)
            pz=round(A[i]['pz'],4)
            theta=round(A[i]['theta'],4)
            phi=round(A[i]['phi'],4)
            f.write(str(E0)+'\t'+str(px)+'\t'+str(py)+'\t'+str(pz)\
                    +'\t'+str(theta)+'\t'+str(phi)+'\n')
            break
        elif (E0<Ecri):
            print(E0)
            E0=round(E0,4)
            px=round(A[i]['px'],4)
            py=round(A[i]['py'],4)
            pz=round(A[i]['pz'],4)
            theta=round(A[i]['theta'],4)
            phi=round(A[i]['phi'],4)
            f.write(str(E0)+'\t'+str(px)+'\t'+str(py)+'\t'+str(pz)\
                    +'\t'+str(theta)+'\t'+str(phi)+'\n')
            i+=1
            continue
        
        theta0=A[i]['theta']
        phi0=A[i]['phi']

        z=gen_z()
        theta1=gen_theta()
        phi1=random.random()*np.pi
    
        E1=(1-z*1/2)*E0 #considering E1 is the more energetic one
        E2=E0-E1
        theta2=-E1*theta1/E2
        phi2=phi1+np.pi

        p1x=E1*np.sin(theta1)*np.cos(phi1)
        p1y=E1*np.sin(theta1)*np.sin(phi1)
        p1z=E1*np.cos(theta1)
        p2x=E2*np.sin(theta2)*np.cos(phi2)
        p2y=E2*np.sin(theta2)*np.sin(phi2)
        p2z=E2*np.cos(theta2)
        J1=R(np.array([E1,p1x,p1y,p1z]),theta0,phi0)
        J2=R(np.array([E2,p2x,p2y,p2z]),theta0,phi0)
        
        theta1=np.arccos(J1[3]/E1)
        phi1=np.arctan(J1[2]/J1[1])
        theta2=np.arccos(J2[3]/E2)
        phi2=np.arctan(J2[2]/J2[1])
                
        B=np.array([(E1,J1[1],J1[2],J1[3],theta1,phi1)],\
                    dtype=partontype)
        C=np.array([(E2,J2[1],J2[2],J2[3],theta2,phi2)],\
                    dtype=partontype)
        A=np.concatenate((A,B,C))
        #print(A)
        i+=1
    f.close()


E1=gen_E0() 
Ecri=20.0 #(Gev)
theta1=random.random()*np.pi
phi1=random.random()*np.pi 
p1x=E1*np.sin(theta1)*np.cos(phi1)
p1y=E1*np.sin(theta1)*np.sin(phi1)
p1z=E1*np.cos(theta1)
A1=np.array([(E1,p1x,p1y,p1z,theta1,phi1)],dtype=partontype)

E2=E1
theta2=np.pi-theta1
phi2=phi1+np.pi
p2x=E2*np.sin(theta2)*np.cos(phi2)
p2y=E2*np.sin(theta2)*np.sin(phi2)
p2z=E2*np.cos(theta2)
A2=np.array([(E2,p2x,p2y,p2z,theta2,phi2)],dtype=partontype)

PartonShower3D(A1,'FinalStates.txt','w')
PartonShower3D(A2,'FinalStates.txt','a')

