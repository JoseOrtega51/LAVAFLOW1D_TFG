# -*- coding: utf-8 -*-
"""
Created on Sun Nov 20 00:08:15 2022

@author: pepeo
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from celluloid import Camera

import sys
EPSILON = sys.float_info.epsilon  # Smallest possible difference.


#Fontsize for figures:
FS1=18
FS2=22
#Linewidth
LWT=2.5


fig=plt.figure(constrained_layout=True)
camera = Camera(fig)
inicio=00
fin=71
x,h,v,zb,T=np.loadtxt("OutputFiles\\datosX_t1.out",unpack=True)
for i in range(inicio,fin):
    x_end,h_end,v_end,zb_end,T_end=np.loadtxt("OutputTransient\\Flow_{:02d}.txt".format(i),unpack=True) 
    plt.plot(x,zb,'k',label="Fondo")
    plt.plot(x,zb+h_end,'b',label="Superficie libre")
    # colors=convert_to_rgb(np.min(T), np.max(T), T_end, colors)
    plt.fill_between(x, zb+h_end, zb, color='blue', alpha=.1)
    plt.plot(x,zb+h,':b',label="Condición inicial")
    plt.xlabel("$x$ (m)")
    plt.ylabel("$z$ (m)")
    plt.xticks()
    plt.yticks()
    camera.snap()
    plt.grid()
anim = camera.animate(blit=True)
anim.save('OutputGIF\\zs.gif')

fig=plt.figure(constrained_layout=True)
camera = Camera(fig)

for i in range(inicio,fin):
    x_end,h_end,v_end,zb_end,T_end=np.loadtxt("OutputTransient\\Flow_{:02d}.txt".format(i),unpack=True)
    plt.plot(x,h_end,'.b',label="Superficie libre")
    # colors=convert_to_rgb(np.min(T), np.max(T), T_end, colors)
    plt.fill_between(x, h_end, 0, color='blue', alpha=.1)
    plt.plot(x,h,':b',label="Condición inicial")
    plt.xlabel("$x$ (m)")
    plt.ylabel("$z$ (m)")
    plt.xticks()
    plt.yticks()
    camera.snap()
    plt.grid()
anim = camera.animate(blit=True)
anim.save('OutputGIF\\h.gif')

fig=plt.figure(constrained_layout=True)
camera = Camera(fig)

for i in range(inicio,fin):
    x_end,h_end,v_end,zb_end,E_end=np.loadtxt("OutputTransient\\Flow_{:02d}.txt".format(i),unpack=True)
    plt.plot(x,v_end,'.b',label="Condición inicial")
    plt.xlabel("$x$ (m)")
    plt.ylabel("$v$ (m/s)")
    plt.xticks()
    plt.yticks()
    plt.ylim(0,1.5)
    camera.snap()
anim = camera.animate(blit=True)
anim.save('OutputGIF\\v.gif')

fig=plt.figure(constrained_layout=True)
camera = Camera(fig)

for i in range(inicio,fin):
    x_end,h_end,v_end,zb_end,T_end=np.loadtxt("OutputTransient\\Flow_{:02d}.txt".format(i),unpack=True)
    plt.plot(x,T_end,'.r')
    plt.xlabel("$x$ (m)")
    plt.ylabel("$T$ (K)")
    plt.xticks()
    plt.yticks()
    camera.snap()
    plt.grid()
    plt.ylim((1000,1500+100))
anim = camera.animate(blit=True)
anim.save('OutputGIF\\T.gif')

