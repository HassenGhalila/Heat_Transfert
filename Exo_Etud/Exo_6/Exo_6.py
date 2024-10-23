#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 18:50:26 2023

@author: hassenghalila

COMMENTS: 
(1) In the plot, recognize that the heat rate for the water is more than 10 times 
    that with oil and 300 times that with air. 
    How do changes in the velocity affect the heat rates for each of the fluids?
    
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


#------------
def Affect():
    global Visc,Pran,Body,Cond,N_Body
    global xmin,Diam,u,u_inf,T_fl,T_s,npt
        
    Re = 5.E+5
    f = open('Data_Exo6.txt', 'r') # lire le fichier

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
    Pres  = float(f.readline().split(':')[0])
    u_inf = float(f.readline().split(':')[0])*1E+0
    u_min = float(f.readline().split(':')[0])*1E+0
    u_max = float(f.readline().split(':')[0])*1E+0
    
    u = np.linspace(u_min,u_max,npt)

    line=f.readline()    
    Visc = []
    Cond = []
    Pran = []
    Body = []
    s = f.readline().split()
    Visc.append(float(s[0])*1.E-6)
    Cond.append(float(s[1]))
    Pran.append(float(s[2]))
    Body.append(s[4])
    line=f.readline()
    N_Body=1
    while True:
        s = line.split()
        if s[0]=='#': break
        Visc.append(float(s[0])*1.E-6)
        Cond.append(float(s[1]))
        Pran.append(float(s[2]))
        Body.append(s[4])
        N_Body += 1
        line=f.readline()
        if not line: break
    f.close()


#------------
def Churchill_Bernstein():
    global Re_D,Nu_D,h_bar,q_dot
    
    Re_D  = np.zeros((N_Body))
    Nu_D  = np.zeros((N_Body))
    h_bar = np.zeros((N_Body))
    q_dot = np.zeros((N_Body))
    for i in range(N_Body):
        Re_D[i] = u_inf*Diam/Visc[i]
        Nu_D[i] = 0.3 + (0.62*np.sqrt(Re_D[i])*Pran[i]**(1/3))/(1+(0.4/Pran[i])**(2/3))**0.25*(1+(Re_D[i]/282000)**(5/8))**(4/5)
        h_bar[i]  = Cond[i]*Nu_D[i]/Diam
        q_dot[i] = h_bar[i]*np.pi*Diam*(T_s - T_fl)

    Re_u = np.zeros((N_Body,npt))
    Nu_u = np.zeros((N_Body,npt))
    h_u  = np.zeros((N_Body,npt))
    q_u  = np.zeros((N_Body,npt))
    for j in range(npt):
        for i in range(N_Body):
            Re_u[i][j] = u[j]*Diam/Visc[i]
            Nu_u[i][j] = 0.3 + (0.62*np.sqrt(Re_u[i][j])*Pran[i]**(1/3))/(1+(0.4/Pran[i])**(2/3))**0.25*(1+(Re_u[i][j]/282000)**(5/8))**(4/5)
            h_u[i][j]  = Cond[i]*Nu_u[i][j]/Diam
            q_u[i][j] = h_u[i][j]*np.pi*Diam*(T_s - T_fl)
            
    return Re_u,Nu_u,h_u,q_u

def Sorties():
    f = open("MyOutPut.txt","w")
    print('\n-----------------------------------------\n')
    print("Results also written in the file : ",f.name)
    print('\n------------------------------------------')

    Output_txt='{0: <10}'.format('Material')+'{0: <13}'.format('  Re_D')+'{0: <10}'.format('Nu_u')+'{0: <14}'.format('h_bar(W/m2/K)')+'{0: <10}'.format('q_dot(W)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=4*'{:6.4e}  '
    for i in range(N_Body):  
        print(f'{Body[i]:10}', 4*'%6.4e   ' %(Re_D[i],Nu_D[i],h_bar[i],q_dot[i]))
        f.write('{:6} : '.format(Body[i]) + MyString.format(Re_D[i],Nu_D[i],h_bar[i],q_dot[i])+'\n')
    f.close()

    PlotFigures()
    
    
#------------
def PlotFigures():    

    Re_u,Nu_u,h_u,q_u = Churchill_Bernstein()
    gs = gridspec.GridSpec(1, 2)    
    fig, CS0 = plt.subplots(figsize=(10, 4))
    CS0 = plt.subplot(gs[0])
    CS1 = plt.subplot(gs[1])
    
    
    
    CS0.plot(u,h_u[0]*100,color='k', alpha=0.8,label='Air (*100)')
    CS0.plot(u,h_u[1],color='g', alpha=0.8,label='Water')
    CS0.plot(u,h_u[2]*10,color='r', alpha=0.8,label='Oil (*10)')
    legend=CS0.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS0.grid()
    CS0.set_title('Churchill Bernstein', fontsize=14,fontstyle='oblique',fontweight="bold")
    CS0.set_xlabel('u (m/s) ')
    CS0.set_ylabel(r'$h (W/m^{2}/K)$')
    
    CS1.plot(u,q_u[0]*100,color='k', alpha=0.8,label='Air (*100)')
    CS1.plot(u,q_u[1],color='g', alpha=0.8,label='Water')
    CS1.plot(u,q_u[2]*10,color='r', alpha=0.8,label='Oil (*10)')
    legend=CS1.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS1.grid()
    CS1.set_title('Churchill Bernstein', fontsize=14,fontstyle='oblique',fontweight="bold")
    CS1.set_xlabel('u (m/s) ')
    CS1.set_ylabel(r'$\dot{q}$ (W)')
    
    plt.tight_layout()    
    plt.show()
    
    

if __name__ == '__main__':
    Affect()
    Churchill_Bernstein()
    Sorties()
