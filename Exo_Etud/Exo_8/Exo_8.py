#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 18:50:26 2023

@author: hassenghalila

COMMENTS:  
(1) Electronic circuit with such measurement sensitivity requires care in its design.
(2) Instruments built on this principle to measure air velocities are called 
    hot-wire anemometers. Generally, the wire diameters are much smaller 
    (3 to 30 µm vs 500 µm of this problem) in order to have faster response times.
(3) What effect would the presence of radiation exchange between the wire 
    and its surroundings have?
    
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

#---------------------
def Affect():
    global Visc,Pran_f,Pran_s,Body_f,Cond_f,N_Body_f
    global Dens_s,Cp,Cond_s,Body_s
    global xmin,Diam_1,Diam_2,C1,m1,C2,m2
    global u,u_inf,T_fl,T_s,T_fin,npt,Correl,Geom
        
    Re = 5.E+5
    f = open('Data_Exo8.txt', 'r') # lire le fichier

    line=f.readline()
    while line[0]=="#":
        line=f.readline()    
    
    npt   = int(line.split(':')[0])
    Sym  = f.readline().split(':')[0]
    xmin  = float(f.readline().split(':')[0])*1E+0
    Diam_1  = float(f.readline().split(':')[0])*1E+0
    Diam_2  = 30.E-6
    s = f.readline().split()
    C1,m1 = float(s[0]),float(s[1])
    s = f.readline().split()
    C2,m2 = float(s[0]),float(s[1])
    Geom = int(f.readline().split(':')[0])
    
    f.readline()
    T_fl     = float(f.readline().split(':')[0])
    T_s     = float(f.readline().split(':')[0])
    T_fin  = float(f.readline().split(':')[0])
    Pres   = float(f.readline().split(':')[0])
    u_inf  = float(f.readline().split(':')[0])*1E+0
    u_min  = float(f.readline().split(':')[0])*1E+0
    u_max  = float(f.readline().split(':')[0])*1E+0
    Correl = int(f.readline().split(':')[0])
    
    u = np.linspace(u_min,u_max,npt)

    line=f.readline()    
    Visc   = []
    Cond_f = []
    Pran_f = []
    Pran_s = []
    Body_f = []
    Dens_s = []
    Cp     = []
    Cond_s = []
    Body_s = []
    s = f.readline().split()
    Visc.append(float(s[0])*1.E-6)
    Cond_f.append(float(s[1]))
    Pran_f.append(float(s[2]))
    Pran_s.append(float(s[3]))
    Body_f.append(s[5])
    line=f.readline()
    N_Body_f =1
    while True:
        s = line.split()
        if s[0][0]=='#': break
        Visc.append(float(s[0])*1.E-6)
        Cond_f.append(float(s[1]))
        Pran_f.append(float(s[2]))
        Pran_s.append(float(s[3]))
        Body_f.append(s[5])
        N_Body_f += 1
        line=f.readline()
        if not line: break
    line=f.readline()
    N_Body_s =0
    while True:
        s = line.split()
        if s[0][0]=='#': break
        Dens_s.append(float(s[0])*1.E-5)
        Cp.append(float(s[1]))
        Cond_s.append(float(s[2]))
        Body_s.append(s[4])
        N_Body_s += 1
        line=f.readline()
        if not line: break
    f.close()  

#---------------------
def Hilpert():
    global Re_D_1,Nu_D_1,h_bar_1,q_dot_1,I_1
    global Re_D_2,Nu_D_2,h_bar_2,q_dot_2,I_2

    #-------- Geometry 1 -------------    
    Re_D_1  = np.zeros((N_Body_f))
    Nu_D_1  = np.zeros((N_Body_f))
    h_bar_1 = np.zeros((N_Body_f))
    q_dot_1 = np.zeros((N_Body_f))
    I_1 = np.zeros((N_Body_f))
    for i in range(N_Body_f):
        Re_D_1[i] = u_inf*Diam_1/Visc[i]
        Nu_D_1[i] = C1*Re_D_1[i]**m1*Pran_f[i]**(1/3)
        h_bar_1[i]  = Cond_f[i]*Nu_D_1[i]/Diam_1
        q_dot_1[i] = h_bar_1[i]*np.pi*Diam_1*(T_s - T_fl)
        I_1[i] = (0.25*np.pi**2*h_bar_1[i]*Diam_1**3/Dens_s[0]*(T_s-T_fl))**0.5
#        print(Re_D_1[i],Nu_D_1[i],h_bar_1[i],q_dot_1[i],I_1[i])

    #-------- Geometry 2 -------------    
    Re_D_2  = np.zeros((N_Body_f))
    Nu_D_2  = np.zeros((N_Body_f))
    h_bar_2 = np.zeros((N_Body_f))
    q_dot_2 = np.zeros((N_Body_f))
    I_2 = np.zeros((N_Body_f))
    for i in range(N_Body_f):
        Re_D_2[i] = u_inf*Diam_2/Visc[i]
        Nu_D_2[i] = C2*Re_D_2[i]**m2*Pran_f[i]**(1/3)
        h_bar_2[i]  = Cond_f[i]*Nu_D_2[i]/Diam_2
        q_dot_2[i] = h_bar_2[i]*np.pi*Diam_2*(T_s - T_fl)
        I_2[i] = (0.25*np.pi**2*h_bar_2[i]*Diam_2**3/Dens_s[0]*(T_s-T_fl))**0.5
#        print(Re_D_2[i],Nu_D_2[i],h_bar_2[i],q_dot_2[i],I_2[i])
        

    Re_u = np.zeros((N_Body_f,npt))
    Nu_u = np.zeros((N_Body_f,npt))
    h_u  = np.zeros((N_Body_f,npt))
    I_u  = np.zeros((N_Body_f,npt))
    if Geom == 0 :
        for j in range(npt):
            for i in range(N_Body_f):
                Re_u[i][j] = u[j]*Diam_1/Visc[i]
                Nu_u[i][j] = C1*Re_u[i][j]**m1*Pran_f[i]**(1/3)
                h_u[i][j]  = Cond_f[i]*Nu_u[i][j]/Diam_1
                I_u[i][j]  = (0.25*np.pi**2*h_u[i][j]*Diam_1**3/Dens_s[0]*(T_s-T_fl))**0.5
    elif Geom == 1 :
        for j in range(npt):
            for i in range(N_Body_f):
                Re_u[i][j] = u[j]*Diam_2/Visc[i]
                Nu_u[i][j] = C2*Re_u[i][j]**m2*Pran_f[i]**(1/3)
                h_u[i][j]  = Cond_f[i]*Nu_u[i][j]/Diam_2
                I_u[i][j]  = (0.25*np.pi**2*h_u[i][j]*Diam_2**3/Dens_s[0]*(T_s-T_fl))**0.5

    return Re_u,Nu_u,h_u,I_u

   
#-----------------------     
def Zhukauskas():
    global Re_D_1,Nu_D_1,h_bar_1,q_dot_1,I_1
    global Re_D_2,Nu_D_2,h_bar_2,q_dot_2,I_2

    #-------- Geometry 1 -------------    
    Re_D_1  = np.zeros((N_Body_f))
    Nu_D_1  = np.zeros((N_Body_f))
    h_bar_1 = np.zeros((N_Body_f))
    q_dot_1 = np.zeros((N_Body_f))
    I_1 = np.zeros((N_Body_f))
    for i in range(N_Body_f):
        Re_D_1[i] = u_inf*Diam_1/Visc[i]
        Nu_D_1[i] = 0.51*Re_D_1[i]**0.5*Pran_f[i]**0.37*(Pran_f[i]/Pran_s[i])**0.25
        h_bar_1[i] = Cond_f[i]*Nu_D_1[i]/Diam_1
        q_dot_1[i] = h_bar_1[i]*np.pi*Diam_1*(T_s - T_fl)
        I_1[i] = (0.25*np.pi**2*h_bar_1[i]*Diam_1**3/Dens_s[0]*(T_s-T_fl))**0.5

    #-------- Geometry 2 -------------    
    Re_D_2  = np.zeros((N_Body_f))
    Nu_D_2  = np.zeros((N_Body_f))
    h_bar_2 = np.zeros((N_Body_f))
    q_dot_2 = np.zeros((N_Body_f))
    I_2 = np.zeros((N_Body_f))
    for i in range(N_Body_f):
        Re_D_2[i] = u_inf*Diam_2/Visc[i]
        Nu_D_2[i] = 0.51*Re_D_2[i]**0.5*Pran_f[i]**0.37*(Pran_f[i]/Pran_s[i])**0.25
        h_bar_2[i]  = Cond_f[i]*Nu_D_2[i]/Diam_2
        q_dot_2[i] = h_bar_2[i]*np.pi*Diam_2*(T_s - T_fl)
        I_2[i] = (0.25*np.pi**2*h_bar_2[i]*Diam_2**3/Dens_s[0]*(T_s-T_fl))**0.5

    Re_u = np.zeros((N_Body_f,npt))
    Nu_u = np.zeros((N_Body_f,npt))
    h_u  = np.zeros((N_Body_f,npt))
    I_u  = np.zeros((N_Body_f,npt))
    if Geom == 0 :
        for j in range(npt):
            for i in range(N_Body_f):
                Re_u[i][j] = u[j]*Diam_1/Visc[i]
                Nu_u[i][j] = 0.51*Re_u[i][j]**0.5*Pran_f[i]**0.37*(Pran_f[i]/Pran_s[i])**0.25
                h_u[i][j]  = Cond_f[i]*Nu_u[i][j]/Diam_1
                I_u[i][j]  = (0.25*np.pi**2*h_u[i][j]*Diam_1**3/Dens_s[0]*(T_s-T_fl))**0.5
    elif Geom == 1 :
        for j in range(npt):
            for i in range(N_Body_f):
                Re_u[i][j] = u[j]*Diam_2/Visc[i]
                Nu_u[i][j] = 0.51*Re_u[i][j]**0.5*Pran_f[i]**0.37*(Pran_f[i]/Pran_s[i])**0.25
                h_u[i][j]  = Cond_f[i]*Nu_u[i][j]/Diam_2
                I_u[i][j]  = (0.25*np.pi**2*h_u[i][j]*Diam_2**3/Dens_s[0]*(T_s-T_fl))**0.5

    return Re_u,Nu_u,h_u,I_u
    

def Sorties():
    f = open("MyOutPut.txt","w")
    print('\n-----------------------------------------------')
    print("Results wrote also in the file : ",f.name)
    print('-------------------------------------------------')

    
    print('----- Results for 1st geometry ----------')     
    f.write('----- Results for 1th geometry ---------- \n')     
    Output_txt='{0: <10}'.format('Material')+'{0: <13}'.format('  Re_D') \
              +'{0: <10}'.format('Nu_u')+'{0: <14}'.format('h_bar(W/m2/K)') \
              +'{0: <10}'.format('q_dot(W)')+'{0: <10}'.format('I(mA)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=4*'{:6.4e}  '+'{:4.2f}  '
    for i in range(N_Body_f):  
        print(f'{Body_f[i]:10}', 4*'%6.4e  ' %(Re_D_1[i],Nu_D_1[i],h_bar_1[i],q_dot_1[i]), '%4.2f  ' %(I_1[i]*1E+3))
        f.write('{:6} : '.format(Body_f[i]) + MyString.format(Re_D_1[i],Nu_D_1[i],h_bar_1[i],q_dot_1[i],I_1[i]*1E+3)+'\n')

    print('----- Results for 2nd geometry ----------')     
    f.write('----- Results for 2nd geometry ---------- \n')     
    Output_txt='{0: <10}'.format('Material')+'{0: <13}'.format('  Re_D') \
              +'{0: <10}'.format('Nu_u')+'{0: <14}'.format('h_bar(W/m2/K)') \
              +'{0: <10}'.format('q_dot(W)')+'{0: <10}'.format('I(mA)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=4*'{:6.4e}  '+'{:4.2f}  '
    for i in range(N_Body_f):  
        print(f'{Body_f[i]:10}', 4*'%6.4e   ' %(Re_D_2[i],Nu_D_2[i],h_bar_2[i],q_dot_2[i]), '%4.2f  ' %(I_2[i]*1E+3))
        f.write('{:6} : '.format(Body_f[i]) + MyString.format(Re_D_2[i],Nu_D_2[i],h_bar_2[i],q_dot_2[i],I_2[i]*1E+3)+'\n')

    print('\n----------------------')
    DI_1 = 0.0025*I_1[i]*1E+6
    DI_1 = ('%5.2f  ' %(0.0025*I_1[i]*1E+6))
    DI_2 = ('%5.2f  ' %(0.0025*I_2[i]*1E+6))
    print('\nTo measure 1% fractional velocity change, a 0.25% fractional change \n'
          'in current must be measured according to DI/I = 0.25DV/V. \n'
          'Then ∆I = 0.0025*I = ' +DI_1+ '(µA) for the first geometry \n'
          'And  ∆I = 0.0025*I = ' +DI_2+ '(µA) for the secund geometry')
    f.close()

    PlotFigures()
    

    
def PlotFigures():    

    if(Correl == 0 ):
        Re_u,Nu_u,h_u,I_u = Hilpert()
        if Geom == 0 :
            Title = r'$Hilpert (D=500 \mu)$'
        elif Geom == 1 :
            Title = r'$Hilpert (D=30 \mu)$'
    elif(Correl == 1):
        Re_u,Nu_u,h_u,I_u = Zhukauskas()
        if Geom == 0 :
            Title = r'$Zhukauskas (D=500 \mu)$'
        elif Geom == 1 :
            Title = r'$Zhukauskas (D=30 \mu)$'

    
    gs = gridspec.GridSpec(1, 2)    
    fig, CS0 = plt.subplots(figsize=(10, 4))
    CS0 = plt.subplot(gs[0])
    CS1 = plt.subplot(gs[1])
    
    
    
    CS0.plot(u,h_u[0],color='k', alpha=0.8,label='Air (T° fluide)')
    CS0.plot(u,h_u[1],color='b', alpha=0.8,label='Air (T° mean)')
    legend=CS0.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS0.grid()
    CS0.set_title(Title, fontsize=14,fontstyle='oblique',fontweight="bold")
    CS0.set_xlabel('u (m/s) ')
    CS0.set_ylabel(r'$h (W/m^{2}/K)$')

    DI_1 = ('%4.2f  ' %((np.max(I_u[0])-np.min(I_u[0]))/np.min(I_u[0])*1E2))    
    DI_2 = ('%4.2f  ' %((np.max(I_u[1])-np.min(I_u[1]))/np.min(I_u[1])*1E2))    
    CS1.plot(u,I_u[0]*1E3,color='k', alpha=0.8,label='DI/I = '+DI_1+'%')
    CS1.plot(u,I_u[1]*1E3,color='b', alpha=0.8,label='DI/I = '+DI_2+'%')
    legend=CS1.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS1.grid()
    CS1.set_title(Title, fontsize=14,fontstyle='oblique',fontweight="bold")
    CS1.set_xlabel('u (m/s) ')
    CS1.set_ylabel('I (mA)')
    
    plt.tight_layout()    
    plt.show()
    
    
if __name__ == '__main__':
    Affect()
    if(Correl == 0):
        Hilpert()
    elif(Correl == 1):
        Zhukauskas()
    Sorties()

