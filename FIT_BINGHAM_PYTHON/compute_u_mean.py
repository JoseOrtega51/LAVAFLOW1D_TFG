# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 22:07:14 2023

@author: pepeo
"""

from lib.Bingham_model import du, tau_y #,mu
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate 
import os
import scipy.optimize as so

T_model=4

#Linear temperature profile given by T hot and T mean
def mu(T,A_mu,B_mu):
    return A_mu*np.exp(B_mu*T)

def T(z,Tbar,h,Th):
    return Th-2./h*(Th-Tbar)*z

def T_lim(z,Tbar,Tf,h):
    DT=(Tbar-Tf)*2./3.
    T=np.zeros(len(z))
    T[z<0.25*h]=Tf+4.*DT/h*z[z<0.25*h]
    T[(z>=0.25*h) & (z<3./4.*h)]=Tf+DT
    T[z>=3./4.*h]=Tf+DT*4*(1-z[z>=3./4.*h]/h)
    return T
    
def T_lim2(z,Tbar,Th,h):
    Tb=4*Tbar-Th*3
    T=np.zeros(len(z))
    T[z<0.25*h]=Tb+(Th-Tb)*z[z<0.25*h]*4/h
    T[(z>=0.25*h) & (z<3./4.*h)]=Th
    T[z>=3./4.*h]=Tb+(Th-Tb)*(h-z[z>=3./4.*h])*4/h
    return T

def T_lim_no_escala(z,Tbar,Tf,h,l):
    Th=(Tbar*h-Tf*l)/(h-l)
    T=np.zeros(len(z))
    T[z<l]=Tf+(Th-Tf)*z[z<l]/l
    T[(z>l) & (z<h-l)]=Th
    T[z>h-l]=Tf+(Th-Tf)*(h-z[z>h-l])/l
    return T
#Folder where to save the data
folder_name="data_u"
if not os.path.exists(folder_name):
    os.mkdir(folder_name)

f=open(folder_name+"\\DATA.txt",'w')
    

#Parameters LAVA NUEVOS
A_tau=0.        #Pa
B_tau=5.6e6    #Pa
C_tau=-0.0058   #K^-1
A_mu=1.77   #Pa s
B_mu=9500.    #K 
T_hot=1500.     #K


# A_tau=0.        #Pa
# B_tau=0.    #Pa
# C_tau=-0.02   #K^-1
# A_mu=0.08      #Pa s
# B_mu=1200.    #K 
# T_hot=600.     #K

# A_tau=0.        #Pa
# B_tau=1000    #Pa
# C_tau=-1./50.   #K^-1
# A_mu=np.exp(-5)      #Pa s
# B_mu=2000.    #K 
# T_hot=323.15     #K
# T_bottom=303.15    #K
# N=20

#Parameters "arbitrarios"
# A_tau=10.        #Pa
# B_tau=2000.    #Pa
# C_tau=-0.02   #K^-1
# A_mu=0.08      #Pa s
# B_mu=1200.    #K 
# tau_b=100.
# T_hot=400.     

# #Benchmark tesis S. MArtinez
# A_tau=0.        #Pa
# B_tau=30431.1    #Pa
# C_tau=-0.00602   #K^-1
# A_mu=3.16    #Pa s
# B_mu=1726.9    #K 
# tau_b=100.
# T_hot=500.  

# #Benchmark D   
# A_tau=0.        #Pa
# B_tau=37000000.    #Pa
# C_tau=-0.01   #K^-1
# A_mu=0.03      #Pa s
# B_mu=9500.    #K 
# tau_b=100.
# T_hot=1050.+273.15     

# #La Palma   
A_tau=0.        #Pa
B_tau=16943587.4    #Pa
C_tau=-0.005   #K^-1
A_mu=0.0555      #Pa s
B_mu=15886.    #K 
T_hot=1500 

#LA PALMA ARTICULO LAY DIFERENTE
A_tau=0.        #Pa
B_tau=2.5e6    #Pa
C_tau=-0.00429   #K^-1
A_mu=2e7      #Pa s
B_mu=-5.02e-3    #K 
T_hot=1500 



T_air=300

N=40
#Parametric sweep
if(T_model==1):
    v_T_mean=np.linspace(T_hot/2,T_hot,N*2)
elif (T_model==2):
    v_T_mean=np.linspace(T_bottom,T_hot,N*2)
elif (T_model==3):
    v_T_mean=np.linspace(T_bottom,T_hot,N*2)
elif (T_model==4):
    v_T_mean=np.loadtxt("out_T.txt")
elif(T_model==5):
    v_T_mean=np.linspace(1/4*T_air+3/4*T_hot,T_hot,N*2)
    
v_tau_b=np.linspace(100,1000000,N)
# v_tau_b=np.linspace(1,100,N*3)
    
#v_h=np.linspace(0.1,10,N)

if(T_model==4):
    v_h=np.array([1.])
    U_mean=np.zeros((len(v_T_mean),1,len(v_tau_b)))
else:
    #v_h=np.linspace(2,10,N)
    v_h=np.array([1.])
    U_mean=np.zeros((N*2,1,len(v_tau_b)))


ind_T=0
ind_h=0
ind_tau=0

for T_mean in v_T_mean:
    for h in v_h:
        for tau_b in v_tau_b:
            #z discretization for trapzoidal integration
            z=np.linspace(0,h,1000)
            
            #Temperture profile evaluation:
            if(T_model==1):
                Temp=T(z,T_mean,h,T_hot)
            elif(T_model==2):
                Temp=T_lim(z,T_mean,T_mean/4,h)
            elif(T_model==3):
                Temp=T_lim_no_escala(z,T_mean,T_bottom,h,1)
            elif(T_model==4):
                z,Temp=np.loadtxt("OutputFiles\\outputData{}.dat".format(np.where(v_T_mean==T_mean)[0][0]),unpack=True)
                Temp=np.flip(Temp[z<1.])
                z=z[z<1.]
            elif(T_model==5):
                Temp=T_lim2(z,T_mean,T_hot,h)
            
            #Velocity profile integration
            u=integrate.cumulative_trapezoid(du(z,Temp,tau_b,h,A_tau,B_tau,C_tau,A_mu,B_mu),z)
            
            
            #Mean velocity
            U_mean[ind_T][ind_h][ind_tau]=integrate.trapezoid(np.concatenate(([0.],u)),z)/h
            np.savetxt(f,np.array([U_mean[ind_T][ind_h][ind_tau],v_T_mean[ind_T],v_h[ind_h],v_tau_b[ind_tau]]).reshape((1,4)))
            
            ind_tau+=1
            
        ind_h+=1
        ind_tau=0
        
    ind_T+=1
    ind_h=0
    ind_tau=0

f.close()    

