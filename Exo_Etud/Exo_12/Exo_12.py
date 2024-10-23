#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 18:50:26 2023

@author: hassenghalila

COMMENTS: 
(1) The lumped capacitance approximation is excellent.
(2) With the surroundings assumed to be at T_sur = T_∞ 
    and a representative emissivity of ε=0.1 for molten aluminum
    the radiation is negligible
    
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#---------------
def Affect():
    global Visc,Pran,Body,Cond,N_Body,Vis_Df,Vis_Ds
    global xmin,Diam,u,u_inf,T_fl,T_s,T_fin,npt
    global Dens_s,Cp,Cond_s,Body_s
        
    Re = 5.E+5
    f = open('Data_Exo12.txt', 'r') # lire le fichier

    line=f.readline()
    while line[0]=="#":
        line=f.readline()    
    
    npt   = int(line.split(':')[0])
    Sym  = f.readline().split(':')[0]
    xmin  = float(f.readline().split(':')[0])*1E+0
    Diam  = float(f.readline().split(':')[0])*1E+0
    f.readline()
    T_fl  = float(f.readline().split(':')[0])
    T_s  = float(f.readline().split(':')[0])
    T_fin  = float(f.readline().split(':')[0])
    Pres  = float(f.readline().split(':')[0])
    u_inf = float(f.readline().split(':')[0])*1E+0
    u_min = float(f.readline().split(':')[0])*1E+0
    u_max = float(f.readline().split(':')[0])*1E+0
    
    u = np.linspace(u_min,u_max,npt)

    line=f.readline()    
    Visc = []
    Cond = []
    Pran = []
    Vis_Df = []
    Vis_Ds = []
    Body   = []
    Dens_s = []
    Cp     = []
    Cond_s = []
    Body_s = []
    s = f.readline().split()
#    print(s)
#    s1,s2 = [float(s[0])*1.E-6, float(s[1])]
#    Dens.append(float(s[0]))
    Visc.append(float(s[0])*1.E-6)
    Cond.append(float(s[1]))
    Pran.append(float(s[2]))
    Vis_Df.append(float(s[3])*1.E-5)
    Vis_Ds.append(float(s[4])*1.E-5)
    Body.append(s[6])
    N_Body=1
    line=f.readline()
    if line[0] !='#':
        while True:
            s = line.split()
            if s[0][0] =='#': break
            Visc.append(float(s[0])*1.E-6)
            Cond.append(float(s[1]))
            Pran.append(float(s[2]))
            Vis_Df.append(float(s[3])*1.E-5)
            Vis_Ds.append(float(s[4])*1.E-5)
            Body.append(s[6])
            N_Body += 1
            line=f.readline()
            if not line: break
#    print(Visc,Cond,Pran,Vis_Df,Vis_Ds,Body,"\n")
    line=f.readline()
    N_Body_s =0
    while True:
        s = line.split()
        if s[0][0]=='#': break
        Dens_s.append(float(s[0]))
        Cp.append(float(s[1]))
        Cond_s.append(float(s[2]))
        Body_s.append(s[4])
        N_Body_s += 1
        line=f.readline()
        if not line: break
#    print(Dens_s,Cp,Cond_s,Body_s,"\n")
    f.close()    

#------------------
def Whitaker():
    global Re_D,Nu_D,h_bar,q_dot,Tof
    
    Re_D  = np.zeros((N_Body))
    Nu_D  = np.zeros((N_Body))
    h_bar = np.zeros((N_Body))
    q_dot = np.zeros((N_Body))
    Tof   = np.zeros((N_Body))
    for i in range(N_Body):
        Re_D[i] = u_inf*Diam/Visc[i]
        Nu_D[i] = 2+(0.4*np.sqrt(Re_D[i])+0.06*Re_D[i]**(2/3))*Pran[i]**0.4*(Vis_Df[i]/Vis_Ds[i])**0.25
        h_bar[i]  = Cond[i]*Nu_D[i]/Diam
        Tof[i] = Dens_s[0]*Cp[0]*Diam/h_bar[i]/6*np.log((T_s-T_fl)/(T_fin-T_fl))
#        q_dot[i] = h_bar[i]*np.pi*Diam**2*(T_s - T_fl)


#------------------
def Sorties():
    f = open("MyOutPut.txt","w")
    print('\n-----------------------------------------')
    print("Results wrote also in the file : ",f.name)

    Output_txt='{0: <10}'.format('Material')+'{0: <16}'.format('  Re_D')+'{0: <10}'.format('Nu_u')+'{0: <15}'.format('h_bar(W/m2/K)')+'{0: <10}'.format('Tof(s)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=4*'{:6.4e}  '
    for i in range(N_Body):  
        print(f'{Body[i]:10}', 4*'%6.4e   ' %(Re_D[i],Nu_D[i],h_bar[i],Tof[i]))
        f.write('{:6} : '.format(Body[i]) + MyString.format(Re_D[i],Nu_D[i],h_bar[i],Tof[i])+'\n')
    f.close()

    print('\n-----------------------------------------')
    print('Maximum thickness of the helium layer needed to ensure the correct T° is :')
    for i in range(N_Body):  
        print(f'{Body[i]:10}',': %6.4f (mm)' %(Tof[i]*u_inf*1E3))

    print('\n-----------------------------------------')
    print('The lumped capacitance approximation and the Bi number :')
    for i in range(N_Body):
        Bi = h_bar[i]*Diam/Cond_s[0]/6
        print(f'{Body[i]:10}',': Bi = %4.2e' %(Bi))
        if Bi  < 0.1 :
            print('Approximation verified for',f'{Body[i]:10}')


if __name__ == '__main__':
    Affect()
    Whitaker()
    Sorties()
