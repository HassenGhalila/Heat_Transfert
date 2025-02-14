#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 18:50:26 2023

@author: hassenghalila

COMMENTS: 
    Note that since Pr>>1, δ>>δt. That is, for the high Prandtl liquids, 
    the velocity boundary layer will be much thicker than the thermal 
	boundary layer
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd


#------------
def Affect():
    global Visc,Pran,Cond,T_fl,x_c,vit
    global xmin,xmax,x,u_inf,npt,ntf,nw,xw,nv
        
    Re = 5.E+5
    f = open('Data_Exo5.txt', 'r') # lire le fichier

    line=f.readline()
    while line[0]=="#":
        line=f.readline()    
    
    npt   = int(line.split(':')[0])
    ntf   = int(f.readline().split(':')[0])
    nw   = int(f.readline().split(':')[0])
    nv   = int(f.readline().split(':')[0])
    Sym  = f.readline().split(':')[0]
    xmin  = float(f.readline().split(':')[0])*1E+0
    xmax  = float(f.readline().split(':')[0])*1E+0
    xw  = float(f.readline().split(':')[0])*1E+0
    f.readline()
    T_fl_min  = float(f.readline().split(':')[0])
    T_fl_max  = float(f.readline().split(':')[0])
    Pres  = float(f.readline().split(':')[0])
    u_inf_min = float(f.readline().split(':')[0])/3.6
    u_inf = float(f.readline().split(':')[0])/3.6
    
    x = np.linspace(xmin,xmax,npt)
    T_fl = np.linspace(T_fl_min,T_fl_max,ntf)
    vit = np.linspace(u_inf_min,u_inf,nv)

    line=f.readline()    
    Visc = []
    Cond = []
    Pran = []
    s = f.readline().split()
    
    MyData_Air = pd.read_excel('TermoPhys.xlsx',sheet_name=0).set_index('T')
    for j in range(ntf):        
        for i in range(len(MyData_Air)-1):
            if T_fl[j] >= MyData_Air.index[i] and T_fl[j] < MyData_Air.index[i+1]:
                a = (MyData_Air.iloc[i+1, 3] - MyData_Air.iloc[i, 3])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 3] - a*MyData_Air.index[i] 
                real_value = a*T_fl[j] + b
                Visc.append(real_value*1.E-6)
                
                a = (MyData_Air.iloc[i+1, 4] - MyData_Air.iloc[i, 4])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 4] - a*MyData_Air.index[i] 
                real_value = a*T_fl[j] + b
                Cond.append(real_value*1.E-3)

                a = (MyData_Air.iloc[i+1, 6] - MyData_Air.iloc[i, 6])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 6] - a*MyData_Air.index[i] 
                real_value = a*T_fl[j] + b
                Pran.append(real_value)

    x_c  = np.zeros((nv,ntf))
    for i in range(ntf):        
        for k in range(nv):        
            x_c[k,i] = Visc[i]*Re/vit[k]
    f.close()    



#-----------------------------------------------
def Thickness():
    global Visc,Pran,Cond,T_fl,x_c,vit
    global xmin,xmax,x,u_inf,npt,ntf,nw,xw,nv
    global h_cv,xxw

    h_cv = np.zeros((nw,nv,ntf))
    xxw  = np.zeros(nw)
    for j in range(ntf):
        for k in range(nv):
            for i in range(nw):
                h_i = 0
                xxw[i] = (i+1)*xw
                if i*xw <= x_c[k,j] and (i+1)*xw < x_c[k,j] :
                    Re_i = vit[k]*i*xw/Visc[j]
                    if i!=0 : h_i = 0.664*Re_i**0.5*Pran[j]**(1/3)*Cond[j]/(i*xw)    
                    Re_f = vit[k]*(i+1)*xw/Visc[j]
                    h_f = 0.664*Re_f**0.5*Pran[j]**(1/3)*Cond[j]/((i+1)*xw) 
                    h_cv[i,k,j] = (h_f*(i+1)*xw - h_i*i*xw)/xw
                else :
                    Re_i = vit[k]*i*xw/Visc[j]
                    if i!=0 : h_i = (0.037*Re_i**(4/5)-871)*Pran[j]**(1/3)*Cond[j]/(i*xw)    
                    Re_f = vit[k]*(i+1)*xw/Visc[j]
                    h_f = (0.037*Re_f**(4/5)-871)*Pran[j]**(1/3)*Cond[j]/((i+1)*xw) 
                    h_cv[i,k,j] = (h_f*(i+1)*xw - h_i*i*xw)/xw



#----------------------------------------------
def Sorties():

    #  Graphics
        
    global h_cv,xxw,nw,nv,ntf,T_fl
    gs = gridspec.GridSpec(2, 2)    
    fig, CS0 = plt.subplots(figsize=(10, 5))
    CS0 = plt.subplot(gs[0])
    CS1 = plt.subplot(gs[1])
    CS2 = plt.subplot(gs[1,:])

    win = 2    
    for i in range(ntf): 
        CS0.plot(vit*3.6,h_cv[win,:,i], alpha=0.8,label='Tf='+'%.2f'%T_fl[i])
    legend=CS0.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS0.grid()
    CS0.set_title('Window n°: %2i'%(win+1), fontsize=14,fontstyle='oblique',fontweight="bold")
    CS0.set_xlabel('v (km/h) ')
    CS0.set_ylabel('h (W/$m^{2}$/K)')

    speed = nv-1
    for i in range(ntf): 
        CS1.plot(xxw,h_cv[:,speed,i], alpha=0.8,label='Tf='+'%.2f'%T_fl[i])
    legend=CS1.legend(loc='upper right', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS1.grid()
    CS1.set_title('U (Km/h) : %.2f'%(vit[speed]*3.6), fontsize=14,fontstyle='oblique',fontweight="bold")
    CS1.set_xlabel('x (m) ')
    CS1.set_ylabel('h (W/$m^{2}$/K)')

    temp = ntf-1
    for i in range(nv): 
        CS2.plot(xxw,h_cv[:,i,temp], alpha=0.8,label='U = '+'%.2f'%(vit[i]*3.6))
    legend=CS2.legend(loc='upper right', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS2.grid()
    CS2.set_title('T film (K): %.2f'%T_fl[temp], fontsize=14,fontstyle='oblique',fontweight="bold")
    CS2.set_xlabel('x (m) ')
    CS2.set_ylabel('h (W/$m^{2}$/K)') 
    
    plt.tight_layout()    
    plt.show()


if __name__ == '__main__':
    Affect()
    Thickness()
    Sorties()

