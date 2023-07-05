# -*- coding: utf-8 -*-
"""
Equations to model the temmperature dependence of the Bingham model parameters and other functions related.
"""

import numpy as np


def tau_y(T,At,Bt,Ct):
    """
    Yield stress dependance with temperature

    Parameters
    ----------
    T : float or array-like
        Temperature [K]
    At : float
        Minimum yield stress [Pa]
    Bt : float
        Decay amplitude [Pa]
    Ct : float
        Decay rate [K^-1]

    Returns
    -------
    array of floats
        Yield stress profile

    """
    return At+Bt*np.exp(Ct*T)

# def mu(T,Au,Bu):
#     """
#     Viscosity dependence with temperature

#     Parameters
#     ----------
#     T : float or array-like
#         Temperature [K]
#     Au : float
#         Infinite temperature limit viscosity [Pa·s]
#     Bu : float
#         Exponential rate [K]

#     Returns
#     -------
#     array of floats
#         Viscosity profile

#     """
#     return Au*np.exp(Bu/T)   

def mu(T,A_mu,B_mu):
    return A_mu*np.exp(B_mu*T)

def du(z,Temp,tau_b,h,At,Bt,Ct,Au,Bu):
    """
    Velocity derivative over the vertical coordinate

    Parameters
    ----------
    z : float or array-like
        Vertical coordinate [m]
    Temp : float or array-like
        Temperature [K]
    tau_b : float or array like
        Bottom stress [Pa]
    h : float
        Fluid height [m]
    At : float
        Minimum yield stress [Pa]
    Bt : float
        Decay amplitude of yield stress [Pa]
    Ct : float
        Decay rate of yield stress[K^-1]
    Au : float
        Infinite temperature limit viscosity [Pa·s]
    Bu : float
        Exponential rate of viscosity [K]

    Returns
    -------
    dudz : array of floats
        Velocity derivative

    """
    tau_umbral=tau_y(Temp,At,Bt,Ct)
    mu_T=mu(Temp,Au,Bu)
    tau=tau_b*(1-z/h)
    dudz=np.zeros(len(z))

    for k in range(0,len(dudz)):
        if(tau[k]>tau_umbral[k]):
            dudz[k]=(tau_b*(1-z[k]/h)-tau_umbral[k])/mu_T[k]
    return dudz 

def u_const_T(z,tau_b,tau_y,mu,h):
    if(tau_b<tau_y):
        U0=0.
    else:
        U0=tau_b*h/2./mu*(1-tau_y/tau_b)**2
    z0=h*(1-tau_y/tau_b)
    U=np.zeros(len(z))
    U[z<=z0]=(tau_b-tau_y)/mu*z[z<=z0]-tau_b/2./mu/h*z[z<=z0]**2
    U[z>z0]=U0
    return U
