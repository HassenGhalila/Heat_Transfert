#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 18:50:26 2023

@author: hassenghalila

COMMENTS: 
(1) The low value of h suggests that heat transfer by free convection may be
    significant and hence that the total loss by convection exceeds 10.3 W.
(2) The surface of the bulb also dissipates heat to the surrounding by radiation. 
    Further, in an actual light bulb, there is also heat loss by conduction through the socket.

"""
import numpy as np

#-------------
def Affect():
    global Visc,Pran,Body,Cond,N_Body,Vis_Df,Vis_Ds
    global xmin,Diam,u,u_inf,T_fl,T_s,npt
        
    Re = 5.E+5
    f = open('Data_Exo10.txt', 'r') # lire le fichier

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
    Vis_Df = []
    Vis_Ds = []
    Body = []
    s = f.readline().split()
    Visc.append(float(s[0])*1.E-6)
    Cond.append(float(s[1]))
    Pran.append(float(s[2]))
    Vis_Df.append(float(s[3])*1.E-7)
    Vis_Ds.append(float(s[4])*1.E-7)
    Body.append(s[6])
    N_Body=1
    line=f.readline()
    if line[0] !='#':
        while True:
            s = line.split()
            if s[0][0]=='#': break
            Visc.append(float(s[0])*1.E-6)
            Cond.append(float(s[1]))
            Pran.append(float(s[2]))
            Body.append(s[4])
            N_Body += 1
            line=f.readline()
            if not line: break
    f.close()    

    
#---------------
def Whitaker():
    global Re_D,Nu_D,h_bar,q_dot
    
    Re_D  = np.zeros((N_Body))
    Nu_D  = np.zeros((N_Body))
    h_bar = np.zeros((N_Body))
    q_dot = np.zeros((N_Body))
    for i in range(N_Body):
        Re_D[i] = u_inf*Diam/Visc[i]
        Nu_D[i] = 2+(0.4*np.sqrt(Re_D[i])+0.06*Re_D[i]**(2/3))*Pran[i]**0.4*(Vis_Df[i]/Vis_Ds[i])**0.25
        h_bar[i]  = Cond[i]*Nu_D[i]/Diam
        q_dot[i] = h_bar[i]*np.pi*Diam**2*(T_s - T_fl)



#------------
def Sorties():
    f = open("MyOutPut.txt","w")
    print('\n-----------------------------------------\n')
    print("Results also written in the file : ",f.name)
    
    print('\n------------------------------------------')

    Output_txt='{0: <10}'.format('Material')+'{0: <15}'.format('  Re_D')+'{0: <11}'.format('Nu_u')+'{0: <14}'.format('h_bar(W/m2/K)')+'{0: <10}'.format('q_dot(W)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=4*'{:6.4e}  '
    for i in range(N_Body):  
        print(f'{Body[i]:10}', 4*'%6.4e   ' %(Re_D[i],Nu_D[i],h_bar[i],q_dot[i]))
        f.write('{:6} : '.format(Body[i]) + MyString.format(Re_D[i],Nu_D[i],h_bar[i],q_dot[i])+'\n')
    f.close()
   
    
if __name__ == '__main__':
    Affect()
    Whitaker()
    Sorties()
