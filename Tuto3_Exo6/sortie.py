#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 21 18:50:26 2018

@author: hassenghalila
"""
import affect as af
import operation as op
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

def ecrire():
    f = open("MyOutPut.txt","w")
    print('\n-----------------------------------------------')
    print("Results wrote also in the file : ",f.name)
    print('-------------------------------------------------')

    
    print('\n----- Results for internal fluide ----------')     
    f.write('\n----- Results for 1th geometry ---------- \n')     
    Output_txt='{0: <10}'.format('Material')+'{0: <13}'.format('  Re_D') \
              +'{0: <10}'.format('Nu_u')+'{0: <14}'.format('h_bar(W/m2/K)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=3*'{:6.4e}  '
    for i in range(af.N_Body_f):  
        print(f'{af.Body_fi[i]:10}', 3*'%6.4e  ' %(op.Re_Di[i],op.Nu_Di[i],op.h_bari[i]))
        f.write('{:6} : '.format(af.Body_fi[i]) + MyString.format(op.Re_Di[i],op.Nu_Di[i],op.h_bari[i])+'\n')

    print('\n----- Results for external fluide ----------')     
    f.write('\n----- Results for 1th geometry ---------- \n')     
    Output_txt='{0: <10}'.format('Material')+'{0: <13}'.format('  Re_D') \
              +'{0: <10}'.format('Nu_u')+'{0: <14}'.format('h_bar(W/m2/K)')
    print(Output_txt)     
    f.write(Output_txt+'\n')
    MyString=3*'{:6.4e}  '
    for i in range(af.N_Body_f):  
        print(f'{af.Body_f[i]:10}', 3*'%6.4e  ' %(op.Re_D[i],op.Nu_D[i],op.h_bar[i]))
        f.write('{:6} : '.format(af.Body_f[i]) + MyString.format(op.Re_D[i],op.Nu_D[i],op.h_bar[i])+'\n')

    MyString=2*'{:6.4e}  '
    Output_txt='  {0: <10}'.format('T_mo')+'  {0: <11}'.format('  T_so')
    print('\n'+Output_txt)     
    f.write('\n'+Output_txt+'\n')
    print(2*'%6.4e  ' %(op.T_mo-273.15,op.T_so-273.15))
    f.write(MyString.format(op.T_mo-273.15,op.T_so-273.15)+'\n')

    f.close()

#    PlotFigures()
    
    
def PlotFigures():    

    if(af.Correl == 0 ):
        Re_u,Nu_u,h_u,I_u = op.Hilpert()
        if af.Geom == 0 :
            Title = r'$Hilpert (D=500 \mu)$'
        elif af.Geom == 1 :
            Title = r'$Hilpert (D=30 \mu)$'
    elif(af.Correl == 1):
        Re_u,Nu_u,h_u,I_u = op.Zhukauskas()
        if af.Geom == 0 :
            Title = r'$Zhukauskas (D=500 \mu)$'
        elif af.Geom == 1 :
            Title = r'$Zhukauskas (D=30 \mu)$'

    
    gs = gridspec.GridSpec(1, 2)    
    fig, CS0 = plt.subplots(figsize=(10, 4))
    CS0 = plt.subplot(gs[0])
    CS1 = plt.subplot(gs[1])
    
    
    
    CS0.plot(af.u,h_u[0],color='k', alpha=0.8,label='Air (T° fluide)')
    CS0.plot(af.u,h_u[1],color='b', alpha=0.8,label='Air (T° mean)')
    legend=CS0.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS0.grid()
    CS0.set_title(Title, fontsize=14,fontstyle='oblique',fontweight="bold")
    CS0.set_xlabel('u (m/s) ')
    CS0.set_ylabel(r'$h (W/m^{2}/K)$')

    DI_1 = ('%4.2f  ' %((np.max(I_u[0])-np.min(I_u[0]))/np.min(I_u[0])*1E2))    
    DI_2 = ('%4.2f  ' %((np.max(I_u[1])-np.min(I_u[1]))/np.min(I_u[1])*1E2))    
    CS1.plot(af.u,I_u[0]*1E3,color='k', alpha=0.8,label='DI/I = '+DI_1+'%')
    CS1.plot(af.u,I_u[1]*1E3,color='b', alpha=0.8,label='DI/I = '+DI_2+'%')
    legend=CS1.legend(loc='upper left', shadow=True, fontsize=8)
    legend.get_frame().set_facecolor('#ffffcc')
    CS1.grid()
    CS1.set_title(Title, fontsize=14,fontstyle='oblique',fontweight="bold")
    CS1.set_xlabel('u (m/s) ')
    CS1.set_ylabel('I (mA)')
    
    plt.tight_layout()    
    plt.show()