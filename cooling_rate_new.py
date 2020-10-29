# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 21:22:16 2020

@author: Vijit Kanjilal
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
#plt.switch_backend('agg')
from scipy import interpolate



t=[0,7,17,33,66,97]
t1=[0.0,1.0,3.0,6.0,12.0,18.0]

dx=1.0
binning=100

kB = 1.3807e-16 #Boltzman's Constant in CGS
mp= 1.67e-24
UNIT_DENSITY= mp

X, Y, Z = 0.7154, 0.2703, 0.0143
mu = 1/(2*X+0.75*Y+0.5625*Z)
mue = 2./(1+X)
mui = 1/(1/mu-1/mue)
gamma = 5./3



x1 = np.zeros(binning, dtype=np.float64)
y1 =  np.zeros(binning, dtype=np.float64)
                
x2 = np.zeros(binning, dtype=np.float64)
y2=  np.zeros(binning, dtype=np.float64)
                
x3 = np.zeros(binning, dtype=np.float64)
y3=  np.zeros(binning, dtype=np.float64)
                
x4 = np.zeros(binning, dtype=np.float64)
y4=  np.zeros(binning, dtype=np.float64)
               
x5 = np.zeros(binning, dtype=np.float64)
y5= np.zeros(binning, dtype=np.float64)

x6 = np.zeros(binning, dtype=np.float64)
y6= np.zeros(binning, dtype=np.float64)

sigma1, sigma2, sigma3, sigma4, sigma5, sigma6 = None, None, None, None, None, None  


cooling = np.loadtxt('cooltable.dat') #solar metallicity
LAMBDA = interpolate.interp1d(cooling[:,0], cooling[:,1], fill_value="extrapolate")

Tmax=-1e20
Tmin=1e20

for i in range(6):
    f = h5py.File('paper_data_new/data.%04d.dbl.h5'%i,'r')       
    T = np.log10(np.array(f['Timestep_%d/vars/T'%t[i]]).flatten())
    Tmax=np.maximum(Tmax,np.max(T))
    Tmin=np.minimum(Tmin,np.min(T))
    #just for the user to know which file is currently in use
    print(i)
    
print(Tmin)
print(Tmax)

Tmax=Tmax+0.01

for i in range(6):
    f = h5py.File('paper_data_new/data.%04d.dbl.h5'%i,'r') 
    rho = (np.array(f['Timestep_%d/vars/rho'%t[i]])).flatten()
    T = np.log10(np.array(f['Timestep_%d/vars/T'%t[i]])).flatten()
    print(i)
    
    #sorting
    T, rho = (list(t) for t in zip(*sorted(zip(T, rho))))
    rho = np.array(rho)
    T = np.array(T)
    #sorting done
    

    bins=np.linspace(Tmin,Tmax,binning+1)
    
    cr = np.zeros(binning, dtype=np.float64)
    
    bins_i=bins[0]
    bins_o=bins[1]

    j=0
    
    for k in range(len(T)):
        ticket=0
        ne = rho[k]*UNIT_DENSITY/(mue*mp)
        ni = rho[k]*UNIT_DENSITY/(mui*mp)
        
        
        while(ticket==0):
           if (bins_i<=T[k]<bins_o):
                cr[j]=cr[j]+(ne*ni*LAMBDA(10**T[k]))
                ticket=1
           else: 
                j=j+1
                bins_i=bins[j]
                bins_j=bins[j+1]
                #ticket=0
            
    
    data_size=len(T)
    print(i)
     
    for j in range(binning): 
        if i==0:
            x1[j]= bins[j]+((bins[j+1]-bins[j])/2)    
            y1[j]= cr[j] 
            sigma1 = np.sqrt(cr/(data_size*(bins[1]-bins[0]) ))
        
        if i==1:
            x2[j]= bins[j]+((bins[j+1]-bins[j])/2)    
            y2[j]= cr[j] 
            sigma2 = np.sqrt(cr/(data_size*(bins[1]-bins[0]) ))
        
        if i==2:
            x3[j]= bins[j]+((bins[j+1]-bins[j])/2)    
            y3[j]= cr[j] 
            sigma3 = np.sqrt(cr/(data_size*(bins[1]-bins[0]) ))
        
        if i==3:
            x4[j]= bins[j]+((bins[j+1]-bins[j])/2)    
            y4[j]= cr[j] 
            sigma4 = np.sqrt(cr/(data_size*(bins[1]-bins[0]) ))
        
        if i==4:
            x5[j]= bins[j]+((bins[j+1]-bins[j])/2)    
            y5[j]= cr[j] 
            sigma5 = np.sqrt(cr/(data_size*(bins[1]-bins[0]) ))
        
        if i==5:
            x6[j]= bins[j]+((bins[j+1]-bins[j])/2)    
            y6[j]= cr[j] 
            sigma6 = np.sqrt(cr/(data_size*(bins[1]-bins[0]) ))
            
fig=plt.figure(figsize=(14,14))

plt.plot(x1,y1,label=r'$t=0$')
plt.fill_between(x1,y1-sigma1,y1+sigma1,alpha=0.4)

plt.plot(x2,y2,label=r'$t\approx t_{cc}$')
plt.fill_between(x2,y2-sigma2,y2+sigma2,alpha=0.4)

plt.plot(x3,y3,label=r'$t\approx 3\, t_{cc}$')
plt.fill_between(x3,y3-sigma3,y3+sigma3,alpha=0.4)

plt.plot(x4,y4,label=r'$t\approx 6\, t_{cc}$')
plt.fill_between(x4,y4-sigma4,y4+sigma4,alpha=0.4)

plt.plot(x5,y5,label=r'$t\approx 12\, t_{cc}$')
plt.fill_between(x5,y5-sigma5,y5+sigma5,alpha=0.4)

plt.plot(x6,y6,label=r'$t\approx 18\, t_{cc}$')
plt.fill_between(x6,y6-sigma6,y6+sigma6,alpha=0.4)

np.savez('data_0_cr.npz',x1,y1,sigma1)
np.savez('data_1_cr.npz',x2,y2,sigma2)
np.savez('data_3_cr.npz',x3,y3,sigma3)
np.savez('data_6_cr.npz',x4,y4,sigma4)
np.savez('data_12_cr.npz',x5,y5,sigma5)
np.savez('data_18_cr.npz',x6,y6,sigma6)


plt.yscale('log')
plt.grid()
plt.legend(loc='best', prop={'size': 24})
plt.xlabel('Temperature',fontsize=22)
plt.ylabel('Cooling Rate',fontsize=22)
plt.tick_params(axis='both', which='major', labelsize=18)
plt.tick_params(axis='both', which='minor', labelsize=14)
plt.tight_layout()  
plt.savefig('Cooling Rate',transparent=True)

