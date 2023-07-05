# -*- coding: utf-8 -*-
"""
Comparison between bingham model with constant temperature and assumed temperature profile

Author: José Ortega Moya
"""

from lib.Bingham_model import tau_y, mu, du, u_const_T
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate 
import os

#Fontsize for figures:
FS1=18
FS2=22
#Linewidth
LWT=2.5

#Folder where to save the plots
folder_name="comparison_Plots"
if not os.path.exists(folder_name):
    os.mkdir(folder_name)
    
T_model=2

#Linear temperature profile given by T hot and T mean
def T(z,Tbar,h,Th):
    return Th-2./h*(Th-Tbar)*z

def T_lim(z,Tbar,Tf,h):
    DT=(Tbar-Tf)*2./3.
    T=np.zeros(len(z))
    T[z<0.25*h]=Tf+4.*DT/h*z[z<0.25*h]
    T[(z>=0.25*h) & (z<3./4.*h)]=Tf+DT
    T[z>=3./4.*h]=Tf+DT*4*(1-z[z>=3./4.*h]/h)
    return T

#Parameters
# A_tau=10.        #Pa
# B_tau=2000.    #Pa
# C_tau=-0.02   #K^-1
# A_mu=0.08      #Pa s
# B_mu=1200.    #K (?)
# tau_b=100.
# T_mean=270.    #K
# T_hot=400.     #K

#Parameters
A_tau=0.        #Pa
B_tau=2.5e6    #Pa
C_tau=-0.006   #K^-1
A_mu=0.1      #Pa s
B_mu=10362.    #K 
T_hot=900.     #K
T_mean=1000.    #K
tau_b=1e6

h=5    #m

# z0=h*(1-tau_y/tau_b)

#z discretization for trapzoidal integration
z=np.linspace(0,h,1000)

#Temperture profile evaluation:
if(T_model==1):
    Temp=T(z,T_mean,h,T_hot)
elif(T_model==2):
    Temp=T_lim(z,T_mean,T_hot,h)

#Velocity profile integration
u=integrate.cumulative_trapezoid(du(z,Temp,tau_b,h,A_tau,B_tau,C_tau,A_mu,B_mu),z)

tau_y_mean=tau_y(T_mean, A_tau, B_tau, C_tau)
mu_mean=mu(T_mean,A_mu,B_mu)
u_T_mean=u_const_T(z, tau_b, tau_y_mean, mu_mean, h)


plt.figure(constrained_layout=True)
plt.plot(u,z[1:len(z)],linewidth=LWT,label=r'$T(z)$')
plt.plot(u_T_mean,z,'--',linewidth=LWT,label=r'$\bar{T}$')
plt.xticks(fontsize=FS1)
plt.yticks(fontsize=FS1)
plt.xlabel(r'$u$ (m/s)',fontsize=FS2)
plt.ylabel(r'$z$ (m)',fontsize=FS2)
plt.legend(fontsize=FS1)
plt.savefig(folder_name +"\\velocity_profile.pdf")
# plt.close()

plt.figure(constrained_layout=True)
plt.plot(Temp,z,linewidth=LWT)
plt.ylabel(r'$z$ (m)',fontsize=FS2)
plt.xlabel(r'$T$ (K)',fontsize=FS2)
plt.xticks(fontsize=FS1)
plt.yticks(fontsize=FS1)
plt.savefig(folder_name +"\\T_profile.pdf")
# plt.close()

tauy=tau_y(Temp,A_tau,B_tau,C_tau)
U_mean=integrate.trapezoid(np.concatenate(([0.],u)),z)/h

plt.figure(constrained_layout=True)
plt.plot(tauy,z,linewidth=LWT,label=r'$\tau_y(T(z))$')
plt.plot(tau_y_mean*np.ones(len(z)),z,'--',linewidth=LWT,label=r'$\tau_y(\bar{T})$')
plt.plot(tau_b*(1-z/h),z,linewidth=LWT,label=r'$\tau$(z)',color='k')
plt.xticks(fontsize=FS1)
plt.yticks(fontsize=FS1)
plt.ylabel(r'$z$ (m)',fontsize=FS2)
plt.xlabel(r'$\tau $ (Pa)',fontsize=FS2)
plt.legend(fontsize=FS1)
plt.savefig(folder_name +"\\stress_profile.pdf")
# plt.close()




