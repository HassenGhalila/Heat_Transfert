#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 18:50:26 2023

@author: hassenghalila

COMMENTS: 
(1) Note that δ≈δt for air, δ>δt for water, δ>>δt for oil, and δ<δt for mercury. 
    As expected, the boundary layer thicknesses increase 
    with increasing distance from the leading edge.
(2) The value of δt for mercury should be viewed as 
    a rough approximation since the expression for δ/δt  
    was derived subject to the approximation that Pr > 0.6
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd

#-------------
def Affect():
    global x,Visc,Pran,Body,N_Body
    global xmin,xmax,u_inf,T_fl,npt
        
    Re = 5.E+5
    f = open('Data_Exo1_V2.txt', 'r') # lire le fichier

    line=f.readline()
    while line[0]=="#":
        line=f.readline()    
    
    npt   = int(line.split(':')[0])
    Sym  = f.readline().split(':')[0]
    xmin  = float(f.readline().split(':')[0])*1E+0
    xmax  = float(f.readline().split(':')[0])*1E+0
    f.readline()
    T_fl  = float(f.readline().split(':')[0])
    Pres  = float(f.readline().split(':')[0])
    u_inf = float(f.readline().split(':')[0])*1E+0
    
    x = np.linspace(xmin,xmax,npt)

    Visc = []
    Pran = []
    Body = []
    Air_OnOff   = int(f.readline().split(':')[0])
    if Air_OnOff == 1 :
        N_Body=1
        MyData_Air = pd.read_excel('TermoPhys.xlsx',sep=';',sheet_name=0).set_index('T')
        Body.append('Air_T_fl')
        Visc.append(MyData_Air.iloc[4, 3]*1.E-6)
        Pran.append(MyData_Air.iloc[4, 6])
        
    Water_OnOff   = int(f.readline().split(':')[0])
    if Water_OnOff == 1 :
        N_Body  += 1
        MyData_Water = pd.read_excel('TermoPhys.xlsx',sep=';',sheet_name=1).set_index('T')
        Body.append('Water_T_fl')
        Visc.append(MyData_Water.iloc[6, 2]*1.E-6/997)
        Pran.append(MyData_Water.iloc[6, 6])
    
    line=f.readline()   
    s = f.readline().split()
    s1,s2 = [float(s[0])*1.E-6, float(s[1])]
    N_Body  += 1
    Visc.append(s1)
    Pran.append(s2)
    Body.append(s[3])
    line=f.readline()
    while True:
        s = line.split()
        s1,s2 = [float(s[0])*1.E-6, float(s[1])]
        N_Body += 1
        Visc.append(s1)
        Pran.append(s2)
        Body.append(s[3])
        line=f.readline()
        if not line: break
    f.close()    


#-------------
def Thickness():
    global Re_l,Thick_Vl,Thick_tl

    Re_x=np.zeros((N_Body,npt))
    Thick_V=np.zeros((N_Body,npt))
    Thick_t=np.zeros((N_Body,npt))
    for j in range(npt):
        for i in range(N_Body):
            Re_x[i][j] = u_inf*x[j]/Visc[i]
            Thick_V[i][j] = 5*x[j]/np.sqrt(Re_x[i][j]+1E-20)
            Thick_t[i][j] = Thick_V[i][j]/Pran[i]**(1/3)

    Re_l=np.zeros((N_Body))
    Thick_Vl=np.zeros((N_Body))
    Thick_tl=np.zeros((N_Body))
    for i in range(N_Body):  
        Re_l[i] = u_inf*xmax/Visc[i]
        Thick_Vl[i] = 5*xmax/np.sqrt(Re_l[i])*1.E+3
        Thick_tl[i] = Thick_Vl[i]/Pran[i]**(1/3)

    return Thick_V,Thick_t


#-------------
def Sorties():
    Thick_V, Thick_t = Thickness()
#    f = open("MyOutPut.txt","w")
#    print('\n-----------------------------------------\n')
#    print("Results also written in the file : ",f.name)
#    MyString=5*'{:7.5f} '
#    for xx,Thick_V0,Thick_V1,Thick_V2,Thick_V3 in zip(x,Thick_V[0][:],Thick_V[1][:],Thick_V[2][:],Thick_V[3][:]):
#        f.write(MyString.format(xx,Thick_V0,Thick_V1,Thick_V2,Thick_V3)+'\n')
#    f.close()
    
    print('\n------------------------------------------')
    print('{0: <10}'.format('Material'),'{0: <10}'.format('  Re_l'),
          '{0: <10}'.format('delta(mm)'),'{0: <10}'.format('delta_t(mm)'))
    for i in range(N_Body):  
        print(f'{Body[i]:10}', 3*'%6.4e ' %(Re_l[i],Thick_Vl[i],Thick_tl[i]))

#   Graphics
        
    Thick_V= Thick_V*1.E+3
    Thick_t= Thick_t*1.E+3
    xx = x*1.E+3
    gs = gridspec.GridSpec(1, 2)    
    fig, CS0 = plt.subplots(figsize=(10, 4))
    CS0 = plt.subplot(gs[0])
    CS1 = plt.subplot(gs[1])
       
    CS0.plot(xx,Thick_V[0],color='k', alpha=0.8,label='Air')
    CS0.plot(xx,Thick_V[1],color='g', alpha=0.8,label='Water')
    CS0.plot(xx,Thick_V[2],color='r', alpha=0.8,label='Oil')
    CS0.plot(xx,Thick_V[3],color='b', alpha=0.8,label='Hg')
    legend=CS0.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS0.grid()
    CS0.set_title('Velocity Boundary Layer', fontsize=14,fontstyle='oblique',fontweight="bold")
    CS0.set_xlabel('x (mm) ')
    CS0.set_ylabel(r'$\delta$ (mm)')
    
    CS1.plot(xx,Thick_t[0],color='k', alpha=0.8,label='Air')
    CS1.plot(xx,Thick_t[1],color='g', alpha=0.8,label='Water')
    CS1.plot(xx,Thick_t[2],color='r', alpha=0.8,label='Oil')
    CS1.plot(xx,Thick_t[3],color='b', alpha=0.8,label='Hg')
    legend=CS1.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS1.grid()
    CS1.set_title('Thermal Boundary Layer', fontsize=14,fontstyle='oblique',fontweight="bold")
    CS1.set_xlabel('x (mm) ')
    CS1.set_ylabel(r'$\delta_{t}$ (mm)')
        
    plt.tight_layout()    
    plt.show()



if __name__ == '__main__':
    Affect()
    Sorties()


