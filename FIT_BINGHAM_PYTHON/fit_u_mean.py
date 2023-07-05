# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:04:56 2023

@author: pepeo
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
styles_path ="..\\"
plt.style.use([styles_path+r"\styles\axis.mplstyle", styles_path+r"\styles\fonts.mplstyle", styles_path+r"\styles\generalSettings.mplstyle"])


pt3d=False

def recta(x,m,n):
    return m*x+n

def pot_desp(x,A,B,n):
    return A+B*x**n

def parabola(x,A,B,C):
    return A*x**3+B*x**2+C

folder_name="data_u"
U,T,h,tau=np.loadtxt(folder_name+"\\DATA.txt",usecols=(0,1,2,3),unpack=True)



#
T_hot=1500.     #K
T_bottom=800.   #K
T_air=300.      #K

N=40

T_model=5



if(T_model==1):
    v_T_mean=np.linspace(T_hot/2,T_hot,N*2)
elif (T_model==2):
    v_T_mean=np.linspace(T_bottom,T_hot,N*2)
elif (T_model==3):
    v_T_mean=np.linspace(T_bottom,T_hot,N*2)
elif (T_model==4):
    v_T_mean=np.loadtxt("out_T.txt")
    v_T_mean=v_T_mean[0:-1]
elif(T_model==5):
    v_T_mean=np.linspace(1/4*T_air+3/4*T_hot,T_hot,N*2)


p1=tau
p2=U/h
p3=T/T_hot

plt.figure(constrained_layout=True)
plt.plot(p3,(p2)/p1,'.')

# Plot 3d
if(pt3d==True):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.scatter(p1,p3,p2)

v_m=np.zeros(len(v_T_mean))
v_n=np.zeros(len(v_T_mean))
ind=0
plt.figure(constrained_layout=True)
for T_cut in v_T_mean:
    coef,cov=so.curve_fit(recta,p1[(T==T_cut) & (p2>0.)],p2[(T==T_cut) & (p2>0.)])
    m,n=coef
    if (ind%3==0):
        plt.plot(p1[T==T_cut],p2[T==T_cut],'.',markersize=10)
       # plt.plot(p1[T==T_cut],recta(p1[T==T_cut],m,n))
    v_m[ind]=m
    v_n[ind]=n
    ind=ind+1
plt.gca().set_prop_cycle(None)
ind=0
for T_cut in v_T_mean:
    if (ind%3==0):
        plt.plot(p1[T==T_cut],recta(p1[T==T_cut],v_m[ind],v_n[ind]))
    ind=ind+1
plt.xlabel("$\\tau_b$ (Pa·s)")
plt.ylabel("$\\overline{U}/h$ (s$^{-1}$)")
plt.savefig("rectas_fit_1.pdf")
# plt.figure()
# plt.plot(np.log(v_T_mean/T_hot),np.log(v_m-0.0016),'.')
# coef,cov=so.curve_fit(recta,np.log(v_T_mean/T_hot),np.log(v_m-0.0016))
# m,n=coef
# plt.plot(np.log(v_T_mean/T_hot),recta(np.log(v_T_mean/T_hot),m,n))
# plt.ylabel("$log(f_1(\Pi_3))$")

plt.figure(constrained_layout=True)
plt.plot(v_T_mean/T_hot,v_m,'.',markersize=10)
param_bounds=([1e-7,-np.inf,-np.inf],[np.inf,np.inf,np.inf])
coef,cov=so.curve_fit(pot_desp,v_T_mean/T_hot,v_m,p0=[1e-6,5e-6,2],bounds=param_bounds)
A1,B1,n1=coef
plt.plot(v_T_mean/T_hot,pot_desp(v_T_mean/T_hot,A1,B1,n1))
plt.ylabel("$f_1(\Pi_3)$")
plt.xlabel("$\\overline{T}/T_h$")
plt.savefig("f1_fit_1.pdf")


plt.figure(constrained_layout=True)
plt.plot(v_T_mean/T_hot,v_n,'.',markersize=10)

# coef,cov=so.curve_fit(pot_desp,v_T_mean/T_hot,v_n,p0=[-5e-2,-0.11,2])
# A2,B2,n2=coef
# plt.plot(v_T_mean/T_hot,pot_desp(v_T_mean/T_hot,A2,B2,n2))
coef,cov=so.curve_fit(parabola,v_T_mean/T_hot,v_n)
A2,B2,n2=coef
plt.plot(v_T_mean/T_hot,parabola(v_T_mean/T_hot,A2,B2,n2))
plt.ylabel("$f_2(\Pi_3)$")
plt.xlabel("$\\overline{T}/T_h$")
plt.savefig("f2_fit_1.pdf")

fig, ax1 = plt.subplots(constrained_layout=True)
ax2 = ax1.twinx()
ax1.plot(v_T_mean/T_hot,v_m*1000,'g.',markersize=10)
ax1.plot(v_T_mean/T_hot,pot_desp(v_T_mean/T_hot,A1,B1,n1)*1000,'g')
ax2.plot(v_T_mean/T_hot,v_n,'b.',markersize=10)
ax2.plot(v_T_mean/T_hot,parabola(v_T_mean/T_hot,A2,B2,n2),'b')
ax1.set_ylabel("$f_1(\\overline{T}) (m/g)$",color='g')
ax2.set_ylabel("$f_2(\\overline{T})$ (s$^{-1}$)",color='b')
ax1.set_xlabel("$\\overline{T}/T_h$")
ax1.grid(False)
ax2.grid(False)
plt.savefig("f1_f2_fit_1.pdf")

epsilon=1e-2
# tau_b_aprox=(p2-pot_desp(T/T_hot,A2,B2,n2))*B_tau/pot_desp(T/T_hot,A1,B1,n1)
tau_b_aprox=(p2-parabola(T/T_hot,A2,B2,n2))/(pot_desp(T/T_hot,A1,B1,n1))


plt.figure(constrained_layout=True)
plt.plot(T[p2>epsilon]/T_hot,tau_b_aprox[p2>epsilon],'s')
plt.plot(T[p2>epsilon]/T_hot,tau[p2>epsilon],'.',markersize=2)
tau_mean=np.mean(tau[p2>epsilon])
SME=np.sum((tau_b_aprox[p2>epsilon]-tau[p2>epsilon])**2)/np.sum((tau_mean-tau[p2>epsilon])**2)
print(np.sqrt(SME))

file1 = open("DATA_FIT_BINGHAM.txt","w")
file1.write("{}\t{}\t{}\t{}\t{}\t{}\t".format(A1,B1,n1,A2,B2,n2))
file1.close()
