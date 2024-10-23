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


#------------
def Affect():
    global Visc,Pran,Body, Dens, Cond, N_Body
    global xmin,xmax,x,u_inf,T_fl,T_s,npt
        
    Re = 5.E+5
    f = open('Data_Exo2.txt', 'r') # lire le fichier

    line=f.readline()
    while line[0]=="#":
        line=f.readline()    
    
    npt   = int(line.split(':')[0])
    Sym  = f.readline().split(':')[0]
    xmin  = float(f.readline().split(':')[0])*1E+0
    xmax  = float(f.readline().split(':')[0])*1E+0
    f.readline()
    T_fl  = float(f.readline().split(':')[0])
    T_s  = float(f.readline().split(':')[0])
    Pres  = float(f.readline().split(':')[0])
    u_inf = float(f.readline().split(':')[0])*1E+0
    
    x = np.linspace(xmin,xmax,npt)

    line=f.readline()    
    Dens = []
    Visc = []
    Cond = []
    Pran = []
    Body = []
    s = f.readline().split()
    Dens.append(float(s[0]))
    Visc.append(float(s[1])*1.E-6)
    Cond.append(float(s[2]))
    Pran.append(float(s[3]))
    Body.append(s[5])
    N_Body = 1
    f.close()    

#------------
def Thickness():
    global Re_l,Thick_Vl,Thick_tl,h_cl,q_dotl,Tau_sl,D_prim

    Re_x=np.zeros((N_Body,npt))
    Thick_V=np.zeros((N_Body,npt))
    Thick_t=np.zeros((N_Body,npt))
    h_c = np.zeros((N_Body,npt))
    q_dot =np.zeros((N_Body,npt))
    Tau_s =np.zeros((N_Body,npt))
    for j in range(npt):
        for i in range(N_Body):
            Re_x[i][j]    = u_inf*x[j]/Visc[i]
            Thick_V[i][j] = 5*x[j]/np.sqrt(Re_x[i][j]+1E-20)
            Thick_t[i][j] = Thick_V[i][j]/Pran[i]**(1/3)
            if(x[j]!=0):
                h_c[i][j]     = Cond[i]/x[j]*0.332*np.sqrt(Re_x[i][j])*Pran[i]**(1/3)
                Tau_s[i][j]   = 0.5*Dens[i]*u_inf**2*0.664/np.sqrt(Re_x[i][j])
            else: 
                h_c[i][j]   = 'Nan'
                Tau_s[i][j] = 'Nan' 
            q_dot[i][j]   = -h_c[i][j]*(T_s - T_f)


    Re_l = np.zeros(N_Body)
    Thick_Vl = np.zeros(N_Body)
    Thick_tl = np.zeros(N_Body)
    h_cl = np.zeros(N_Body)
    q_dotl = np.zeros(N_Body)
    Tau_sl = np.zeros(N_Body)
    D_prim = np.zeros(N_Body)
    # Reynolds number, and Boundary layers
    for i in range(N_Body): 
        Re_l[i] = u_inf*xmax/Visc[i]
        Thick_Vl[i] = 5*xmax/np.sqrt(Re_l[i])*1.E+3
        Thick_tl[i] = Thick_Vl[i]/Pran[i]**(1/3)
    # The local convection coefficient, Eq. 7.23, and heat flux at x = L 
    for i in range(N_Body): 
        h_cl[i]   = Cond[i]/xmax*0.332*np.sqrt(Re_l[i])*Pran[i]**(1/3)
        q_dotl[i] = -h_cl[i]*(T_s - T_f)
    # The local shear stress (Eq. 7.20) and the drag force per unit width
    for i in range(N_Body):
        Tau_sl[i] = 0.5*Dens[i]*u_inf**2*0.664/np.sqrt(Re_l[i])
        D_prim[i] = 2*xmax*Tau_sl[i]

    return Thick_V,Thick_t,h_c,q_dot,Tau_s


#------------
def Sorties():
    Thick_V, Thick_t,h_c,q_dot,Tau_s = Thickness()
    f = open("MyOutPut.txt","w")
    print('\n-----------------------------------------\n')
    print("Results also written in the file : ",f.name)
    MyString=6*'{:7.5f} '
    for xx,Thick_V0,Thick_t0,h_c0,q_dot0,Tau_s0 in zip(x,Thick_V[0][:],Thick_t[0][:],h_c[0][:],q_dot[0][:],Tau_s[0][:]):
        f.write(MyString.format(xx,Thick_V0,Thick_t0,h_c0,q_dot0,Tau_s0)+'\n')
    f.close()
    
    print('\n------------------------------------------')
    print('{0: <10}'.format('Material'),'{0: <10}'.format('  Re_l'),
          '{0: <10}'.format('delta(mm)'),'{0: <10}'.format('delta_t(mm)'))
    for i in range(N_Body):  
        print(f'{Body[i]:10}', 3*'%6.4e ' %(Re_l[i],Thick_Vl[i],Thick_tl[i]))
    print('\n------------------------------------------')
    for i in range(N_Body):
        print('{0: <10}'.format('Material'),'{0: <10}'.format('  h_c_l'),
              '{0: <10}'.format('q_dot_l'))
        print(f'{Body[i]:10}', 2*'%6.4e ' %(h_cl[i],q_dotl[i]))
    print('\n------------------------------------------')
    for i in range(N_Body):
        print('{0: <10}'.format('Material'),'{0: <10}'.format('  Tau_sl'),
              '{0: <10}'.format('D_prim'))
        print(f'{Body[i]:10}', 2*'%6.4e ' %(Tau_sl[i],D_prim[i]))

#  Graphics
        
    Thick_V, Thick_t, h_c, q_dot, Tau_s = Thickness()
    Thick_V= Thick_V*1.E+0
    Thick_t= Thick_t*1.E+0
    gs = gridspec.GridSpec(2, 2)    
    fig, CS0 = plt.subplots(figsize=(10, 8))
    CS0 = plt.subplot(gs[0])
    CS1 = plt.subplot(gs[1])
    CS2 = plt.subplot(gs[2])
       
    CS0.plot(x,Thick_V[0],color='k', alpha=0.8,label='Oil')
    legend=CS0.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS0.grid()
    CS0.set_title('Velocity Boundary Layer', fontsize=14,fontstyle='oblique',fontweight="bold")
    CS0.set_xlabel('x (m) ')
    CS0.set_ylabel(r'$\delta$ (m)')
    
    CS1.plot(x,Thick_t[0],color='k', alpha=0.8,label='Oil')
    legend=CS1.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS1.grid()
    CS1.set_title('Thermal Boundary Layer', fontsize=14,fontstyle='oblique',fontweight="bold")
    CS1.set_xlabel('x (m) ')
    CS1.set_ylabel(r'$\delta_{t}$ (m)')
    
    CS2.plot(x,Thick_V[0]*1.E+4,color='k', alpha=0.8,label='HBL')
    CS2.plot(x,h_c[0]*1E+2,color='g', alpha=0.8,label=r'$h_{c}$')
    CS2.plot(x,q_dot[0],color='r', alpha=0.8,label='Heat rate ($q_{dot}$)')
    CS2.plot(x,Tau_s[0]*1.E+4,color='b', alpha=0.8,label=r'$\tau_{s}$')
    legend=CS2.legend(loc='upper right', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS2.grid()
    CS2.set_title('Thermal parameters', fontsize=14,fontstyle='oblique',fontweight="bold")
    CS2.set_xlabel('x (m) ')
    CS2.set_ylabel(r'$\delta_{t}*1E4$, $h_{c}*1E2$, $q_{dot}$')
    CS2.set_ylim(0,5000)
   
    plt.tight_layout()    
    plt.show()


if __name__ == '__main__':
    Affect()
    Sorties()

