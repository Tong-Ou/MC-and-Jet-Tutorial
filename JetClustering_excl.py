# -*- coding: utf-8 -*-
"""
Created on Fri May 25 21:35:18 2018

@author: outong
"""

import numpy as np
import numpy.random as random
import time

partontype=np.dtype({'names':['energy','px','py','pz','theta','phi'],\
                     'formats':['f','f','f','f','f','f']})

Ecri=20 # (Gev) Energy threshold 
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
def PartonShower3D(A):
    i=0
    FS=[]
    while 1:
        E0=A[i]['energy']
        if (i==len(A)-1)and(E0<Ecri):
           # print(E0)
            FS.append(A[i])
            break
        elif (E0<Ecri):
           # print(E0)
            FS.append(A[i])
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
#    print(FS)
    return FS

'''
Generate two partons initially coming out 
from the pp collision.
'''
def start():
    E1=gen_E0() 
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
    
    return A1,A2


'''
Parton Shower Generation
-----------------------------------------------------------
-----------------------------------------------------------
Jet Clustering
'''

def pT(px,py):
    return np.sqrt(px**2+py**2)

def Delta(theta1,phi1,theta2,phi2):
    eta1=-np.log(abs(np.tan(theta1/2)))
    eta2=-np.log(abs(np.tan(theta2/2)))
    return np.sqrt((phi1-phi2)**2+(eta1-eta2)**2)

def JetCluster_excl(p,R,A,dcut):
    A0=A
    jet=[] 
    pjet2=[]
    n=len(A) 
    pjet=np.zeros([n,2],dtype=partontype) 
    
    for i in range(n):
        pjet[i,0]=A[i]
        
    diB=np.zeros(n)
    dij=np.zeros([n,n])
    for i in range(n):
        pxi=A[i]['px']
        pyi=A[i]['py']
        thetai=A[i]['theta']
        phii=A[i]['phi']
        
        pti=pT(pxi,pyi)
        diB[i]=pti**(2*p)
        if i==n-1:
            break
        for j in range(i+1,n):
            pxj=A[j]['px']
            pyj=A[j]['py']
            thetaj=A[j]['theta']
            phij=A[j]['phi']
            
            ptj=pT(pxj,pyj)
            pt=min(pti**(2*p),ptj**(2*p))
            delta=Delta(thetai,phii,thetaj,phij)
            dij[i,j]=pt*(delta**2)/(R**2)
            
    while len(A)>1:
        n=len(A)
        
        # To find the minimum among diB and dij
        t1=diB[0]
        m1=0
        for i in range(n):
            if diB[i]<t1:
               t1=diB[i]
               m1=i
        
        t2=dij[0,1]
        m2=[0,1]
        for i in range(n):
            if i==n-1:
                break
            for j in range(i+1,n):
                if dij[i,j]<t2:
                    t2=dij[i,j]
                    m2=[i,j]
                    #m2[0]<m2[1]

        #clustering is stopped when all dij and diB are above dcut
        if min(t1,t2)>dcut:
            break
                  
        if t1<t2:
            if not(A[m1] in A0):
                jet.append(A[m1])
                pjet2.append(pjet[m1])
            A=np.concatenate((A[:m1],A[m1+1:]))
            pjet=np.delete(pjet,m1,0)
            diB=np.concatenate((diB[:m1],diB[m1+1:]))
            dij=np.delete(dij,m1,0)
            dij=np.delete(dij,m1,1)
        else:
            r=m2[0]
            s=m2[1]
            
            #put the two highest pT constituents of a jet into pjet
            px1=pjet[r,0]['px']
            py1=pjet[r,0]['py']
            px2=pjet[r,1]['px']
            py2=pjet[r,1]['py']
            pt1=pT(px1,py1)
            pt2=pT(px2,py2)
            
            px3=pjet[s,0]['px']
            py3=pjet[s,0]['py']
            px4=pjet[s,1]['px']
            py4=pjet[s,1]['py']
            pt3=pT(px3,py3)
            pt4=pT(px4,py4)
            if min(pt3,pt4)>max(pt1,pt2):
                pjet[r,0]=pjet[s,0]
                pjet[r,1]=pjet[s,1]
            elif max(pt3,pt4)>min(pt1,pt2):
                j1=int(pt1>pt2)
                j2=int(pt3<pt4)
                pjet[r,j1]=pjet[s,j2]

            pjet=np.delete(pjet,s,0)

           #calculate the four momenta of the combined pseudojet
            A[r]['energy']+=A[s]['energy']
            A[r]['px']+=A[s]['px']
            A[r]['py']+=A[s]['py']
            A[r]['pz']+=A[s]['pz']
            A[r]['theta']=np.arccos(A[r]['pz']/A[r]['energy'])
            A[r]['phi']=np.arctan(A[r]['py']/A[r]['px'])
            
            A=np.concatenate((A[:s],A[s+1:]))
            
            #recalculate the diB and dij matrices
            diB[r]=pT(A[r]['px'],A[r]['py'])**(2*p)
            diB=np.concatenate((diB[:s],diB[s+1:]))
            dij=np.delete(dij,s,0)
            dij=np.delete(dij,s,1)
            
            ptr=pT(A[r]['px'],A[r]['py'])
            thetar=A[r]['theta']
            phir=A[r]['phi']
            for i in range(r):
                pti=pT(A[i]['px'],A[i]['py'])
                thetai=A[i]['theta']
                phii=A[i]['phi']
                pt=min(pti**(2*p),ptr**(2*p))
                delta=Delta(thetai,phii,thetar,phir)
                dij[i,r]=pt*(delta**2)/(R**2)
            if r<len(A)-1:
                for j in range(r+1,len(A)):
                    ptj=pT(A[j]['px'],A[j]['py'])
                    thetaj=A[j]['theta']
                    phij=A[j]['phi']
                    pt=min(ptj**(2*p),ptr**(2*p))
                    delta=Delta(thetaj,phij,thetar,phir)
                    dij[r,j]=pt*(delta**2)/(R**2)
                             
    njet=len(jet)
    
    '''
    #find the highest energy Jet to calculate its pseudomass
    ee=jet[0]['energy']
    mm=0
    for i in range(njet):
       if  jet[i]['energy']>ee:
           ee=jet[i]['energy']
           mm=i
    
    #calculate the pseudomass
    E1=pjet2[mm][0]['energy']
    E2=pjet2[mm][1]['energy']
    theta1=pjet2[mm][0]['theta']
    phi1=pjet2[mm][0]['phi']
    theta2=pjet2[mm][1]['theta']
    phi2=pjet2[mm][1]['phi']
    if theta2==0:
        
        pmass=0
    else:
        pmass=E1*E2*Delta(theta1,phi1,theta2,phi2)    
    '''    
    
    return njet


A1,A2=start()
#A1=np.array([(579.34320068,-455.52883911,350.43041992,\
 #             -73.01068115,1.69715548,2.4858644)],dtype=partontype)
#A2=np.array([(579.34320068,455.52883911,-350.43041992,\
  #            73.01068115,1.44443715,5.62745714)],dtype=partontype)
#f0=open('A1A2.txt','w')
#f0.write(str(A1)+'\n'+str(A2))
#f0.close()

be=time.clock()
f=open('njet_p0_R1.0_excl.txt','w')
for i in range(1000):
    FS1=PartonShower3D(A1)
    FS2=PartonShower3D(A2)
    FS=np.concatenate((FS1,FS2))
    nj=JetCluster_excl(0,1.0,FS,5000)
    f.write(str(nj)+'\n')
    if i%100==0:
        print(i)
f.close()

end=time.clock()
print('time=',end-be)    
