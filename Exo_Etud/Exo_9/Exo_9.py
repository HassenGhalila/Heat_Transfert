#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 18:50:26 2023

@author: hassenghalila

COMMENTS:  
On Biot number Bi = 0,056 << 0.1 => the lumped capacitance method can be applied
- Measures the fraction of conductive and convective flux Conv/Cond
- It's a dimensionless number used in transient heat transfer calculations.
- If it is small in front of 1 (we will often use Bi < 0.1), 
  this means that the internal resistance is negligible, and therefore that 
  the temperature can be considered uniform inside the body.
- It should not be confused with the Nusselt number, which employs the 
  thermal conductivity of the fluid rather than that of the body

"""

import numpy as np

#-------------
def Affect():
    global Visc,Pran,Body_f,Cond_f,N_Body_f
    global Dens_s,Cp,Cond_s,Body_s
    global xmin,Diam_1,Diam_2,C1,m1,C2,m2
    global u,u_inf,T_fl,T_s,T_fin,npt
        
    Re = 5.E+5
    f = open('Data_Exo9.txt', 'r') # lire le fichier

    line=f.readline()
    while line[0]=="#":
        line=f.readline()    
    
    npt   = int(line.split(':')[0])
    Sym  = f.readline().split(':')[0]
    xmin  = float(f.readline().split(':')[0])*1E+0
    Diam_1  = float(f.readline().split(':')[0])*1E+0
    Diam_2  = Diam_1*np.sqrt(2)
    s = f.readline().split()
    C1,m1 = float(s[0]),float(s[1])
    s = f.readline().split()
    C2,m2 = float(s[0]),float(s[1])
    
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
    Visc   = []
    Cond_f = []
    Pran   = []
    Body_f = []
    Dens_s = []
    Cp     = []
    Cond_s = []
    Body_s = []
    s = f.readline().split()
    Visc.append(float(s[0])*1.E-6)
    Cond_f.append(float(s[1]))
    Pran.append(float(s[2]))
    Body_f.append(s[4])
    line=f.readline()
    N_Body_f =1
    while True:
        s = line.split()
        if s[0][0]=='#': break
        Visc.append(float(s[0])*1.E-6)
        Cond_f.append(float(s[1]))
        Pran.append(float(s[2]))
        Body_f.append(s[4])
        N_Body_f += 1
        line=f.readline()
        if not line: break
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
    f.close()  

	
#-------------
def Hilpert():
    global Re_D_1,Nu_D_1,h_bar_1,q_dot_1,Time_1
    global Re_D_2,Nu_D_2,h_bar_2,q_dot_2,Time_2

    #-------- Geometry 1 -------------    
    Re_D_1  = np.zeros((N_Body_f))
    Nu_D_1  = np.zeros((N_Body_f))
    h_bar_1 = np.zeros((N_Body_f))
    q_dot_1 = np.zeros((N_Body_f))
    Time_1 = np.zeros((N_Body_f))
    for i in range(N_Body_f):
        Re_D_1[i] = u_inf*Diam_1/Visc[i]
        Nu_D_1[i] = C1*Re_D_1[i]**m1*Pran[i]**(1/3)
        h_bar_1[i]  = Cond_f[i]*Nu_D_1[i]/Diam_1
        q_dot_1[i] = h_bar_1[i]*np.pi*Diam_1*(T_s - T_fl)
        Time_1[i] = 0.25*Dens_s[0]*Cp[0]*Diam_1*np.log((T_s-T_fl)/(T_fin-T_fl))/3600

    #-------- Geometry 2 -------------    
    Re_D_2  = np.zeros((N_Body_f))
    Nu_D_2  = np.zeros((N_Body_f))
    h_bar_2 = np.zeros((N_Body_f))
    q_dot_2 = np.zeros((N_Body_f))
    Time_2 = np.zeros((N_Body_f))
    for i in range(N_Body_f):
        Re_D_2[i] = u_inf*Diam_2/Visc[i]
        Nu_D_2[i] = C2*Re_D_2[i]**m2*Pran[i]**(1/3)
        h_bar_2[i]  = Cond_f[i]*Nu_D_2[i]/Diam_2
        q_dot_2[i] = h_bar_2[i]*np.pi*Diam_2*(T_s - T_fl)
        Time_2[i] = 0.25*Dens_s[0]*Cp[0]*Diam_2*np.log((T_s-T_fl)/(T_fin-T_fl))/3600

		
#-------------
def Sorties():
    f = open("MyOutPut.txt","w")
    print('\n-----------------------------------------\n')
    print("Results also written in the file : ",f.name)
    print('\n-----------------------------------------\n')

    
    print('----- Results for 1th geometry ----------')     
    f.write('----- Results for 1th geometry ---------- \n')     
    Output_txt='{0: <10}'.format('Material')+'{0: <13}'.format('  Re_D') \
              +'{0: <10}'.format('Nu_u')+'{0: <14}'.format('h_bar(W/m2/K)') \
              +'{0: <10}'.format('q_dot(W)')+'{0: <10}'.format('Time(H)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=4*'{:6.4e}  '+'{:4.2f}  '
    for i in range(N_Body_f):  
        print(f'{Body_f[i]:10}', 4*'%6.4e  ' %(Re_D_1[i],Nu_D_1[i],h_bar_1[i],q_dot_1[i]), '%4.2f  ' %(Time_1[i]))
        f.write('{:6} : '.format(Body_f[i]) + MyString.format(Re_D_1[i],Nu_D_1[i],h_bar_1[i],q_dot_1[i],Time_1[i])+'\n')

    print('\n----------------------\n')
    print('----- Results for 2nd geometry ----------')     
    f.write('----- Results for 2nd geometry ---------- \n')     
    Output_txt='{0: <10}'.format('Material')+'{0: <13}'.format('  Re_D') \
              +'{0: <10}'.format('Nu_u')+'{0: <14}'.format('h_bar(W/m2/K)') \
              +'{0: <10}'.format('q_dot(W)')+'{0: <10}'.format('Time(H)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=4*'{:6.4e}  '+'{:4.2f}  '
    for i in range(N_Body_f):  
        print(f'{Body_f[i]:10}', 4*'%6.4e   ' %(Re_D_2[i],Nu_D_2[i],h_bar_2[i],q_dot_2[i]), '%4.2f  ' %(Time_2[i]))
        f.write('{:6} : '.format(Body_f[i]) + MyString.format(Re_D_2[i],Nu_D_2[i],h_bar_2[i],q_dot_2[i],Time_2[i])+'\n')

    f.close()


if __name__ == '__main__':
    Affect()
    Hilpert()
    Sorties()
