# -*- coding: utf-8 -*-
"""
Created on Sat Jul  1 19:30:42 2023

@author: pepeo
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so

def exponential(T,a,n):
    return a*(1501-T)**n

styles_path = "..\\"
plt.style.use([styles_path+r"\styles\axis.mplstyle", styles_path+r"\styles\fonts.mplstyle", styles_path+r"\styles\generalSettings.mplstyle"])



T_mean=np.loadtxt("out_T.txt")
N=len(T_mean)
T_s=np.zeros(N)

for i in range(N):
    x,T=np.loadtxt("OutputFiles\\outputData{}.dat".format(i),unpack=True)
    T_s[i]=T[0]

plt.figure(constrained_layout=True)

coef,cov=so.curve_fit(exponential,T_mean[0:N],T_s)
A,B=coef
plt.plot(T_mean[0:N],T_s,'o',markersize=7)
plt.plot(np.linspace(400,1500,1000),exponential(np.linspace(400,1500,1000),A,B),'--',color='k',linewidth=2.5)

plt.xlabel("$\overline{T}$ (K)")
plt.ylabel("$T_s$ (K)")
plt.savefig("Ts_fit.pdf")