# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 23:29:52 2020

@author: Vijit Kanjilal
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate

CHI = np.logspace(1.5, 2.5, 100)
M = np.linspace(0.5, 1.5, 100)
cooling = np.loadtxt('cooltable.dat') #solar metallicity
LAMBDA = interpolate.interp1d(cooling[:,0], cooling[:,1])


X1, X2 = np.meshgrid(CHI, M)

#Constants
kB = 1.3807e-16 #Boltzman's Constant in CGS
mp = 1.6726231e-24 #Mass of a Proton in CGS
GAMMA = 5./3 #Specific Heat Ratio for an Ideal Gas

#Problem Constants
mu = 0.672442
Tcl = 1.e4 #K
ncl = 0.1 # particles per cm^3
T_hot = CHI*Tcl
LAMBDA_HOT= LAMBDA(T_hot) #erg cm3 s-1    #LAMBDA at T_hot #GET IT FROM COOLTABLE.DAT
Tmix= np.sqrt(Tcl*T_hot) #K
LAMBDA_MIX = LAMBDA(Tmix) #erg cm3 s-1    #LAMBDA at T_mix #GET IT FROM COOLTABLE.DAT
ALPHA = 1
n_hot=ncl/CHI


#Normalized Quantities
Tcl_4 = Tcl/1e4 #K
P3 = (ncl*Tcl)/1e3 #cm-3 K 
CHI_100=(X1)/100
LAMBDA_HOT_N21_4 = LAMBDA_HOT/(10**-21.4) #erg cm3 s-1
LAMBDA_MIX_N21_4 = LAMBDA_MIX/(10**-21.4) #erg cm3 s-1


cs = np.sqrt(GAMMA*kB/(mu*mp))*np.sqrt(X1*Tcl)/1e5 #kms^-1
vCl=X2*cs
R_Li = (15.4*(Tcl_4**(12/13)))*(X2**(4/13))*((CHI_100)**(20/13))/((ncl/0.1)*(LAMBDA_HOT_N21_4**(10/13)))
 
#print(R_Li)

Rgo = (2 * (Tcl_4**(5/2)) * X2 * CHI_100 )/(P3*LAMBDA_MIX_N21_4*ALPHA)
#print(Rgo)

#Plotting Area

#Gronke Oh et. al.
fig = plt.figure(figsize=(20,20))
CS10 = plt.contour(X1,X2,Rgo,13,colors='w')
CS11 = plt.contourf(X1,X2,Rgo,50)
cb=plt.colorbar()
cb.set_label(label='Critical Radius (pc)',size=30)
cb.ax.tick_params(labelsize=30)
plt.clabel(CS10,inline=10,fontsize=32,colors='w')
plt.xlabel(r'$\chi$',fontsize=70)
plt.ylabel('Mach No.',fontsize=65)
#plt.title('Gronke Oh Cloud Radius w.r.t. CHI and Mach No.',fontsize=40)
plt.tick_params(axis='both', which='major', labelsize=50, direction="out", pad=15)
plt.tick_params(axis='both', which='minor', labelsize=48, direction="out", pad=15)
plt.grid()
plt.savefig('Gronke_Oh.pdf',transparent =True, bbox_inches='tight')
plt.close()

#Li et. al.
fig = plt.figure(figsize=(20,20))
CS10 = plt.contour(X1,X2,R_Li,10,colors='w')
CS11 = plt.contourf(X1,X2,R_Li,50)
cb=plt.colorbar()
cb.ax.tick_params(labelsize=30)
cb.set_label(label='Critical Radius (pc)',size=30)
plt.clabel(CS10,inline=10,fontsize=32,colors='w')
plt.xlabel(r'$\chi$',fontsize=70)
plt.ylabel('Mach No.',fontsize=65)
#plt.title('Li Cloud Radius w.r.t. CHI and Mach No.',fontsize=40)
plt.tick_params(axis='both', which='major', labelsize=50, direction="out", pad=15)
plt.tick_params(axis='both', which='minor', labelsize=48, direction="out", pad=15)
plt.grid()
plt.savefig('Li.pdf',transparent =True, bbox_inches='tight')
plt.close()

#Ratio 
Z = np.log10(Rgo/R_Li)
fig = plt.figure(figsize=(20,20))
CS10 = plt.contour(X1, X2, Z, 8, colors='w')
CS11 = plt.contourf(X1, X2, Z, 100)
cb=plt.colorbar()
cb.ax.tick_params(labelsize=30)
cb.set_label(label=r'$log\frac{R_{GO}}{R_{Li}}$',size=40)
plt.clabel(CS10, inline=1, fontsize=32, colors='w')
#plt.plot([100, 50, 100, 300, 300, 50], [1.5, 1.5, 2.0, 1.5, 2.0, 2.0], 'o', color='red')
plt.xlabel(r'$\chi$',fontsize=70)
plt.ylabel('Mach No.',fontsize=65)
#plt.xscale('log')
#plt.title(r'$log\frac{R_{GO}}{R_{LH}}$', fontsize=20, y=1.03)
#plt.show()
plt.tick_params(axis='both', which='major', labelsize=50, direction="out", pad=15)
plt.tick_params(axis='both', which='minor', labelsize=48, direction="out", pad=15)
plt.grid()
plt.savefig('Ratio.pdf',transparent =True, bbox_inches='tight')
plt.close()