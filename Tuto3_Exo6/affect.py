#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Friday Oct 13 13:55:20 2023

@author: hassenghalila
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
import sys
import os

#-------------
def Affect():
    global x,Visc,Pran,Body,Cond,Cp,Beta,Alfa,N_Body
    global xmin,xmax,Power,Diam,u_inf,T_fl,T_si,T_fi,npt
    global Ra, RayD, Air_OnOff, Water_OnOff
        
    Ra = 1.E+12
    f = open('Data_Exo6.txt', 'r') # lire le fichier

    line=f.readline()
    while line[0]=="#":
        line=f.readline()    
    
    npt   = int(line.split(':')[0])
    Sym  = f.readline().split(':')[0]
    xmin  = float(f.readline().split(':')[0])*1E+0
    xmax  = float(f.readline().split(':')[0])*1E+0
    Diam  = float(f.readline().split(':')[0])*1E+0
    f.readline()
    T_fl  = float(f.readline().split(':')[0])
    T_si  = float(f.readline().split(':')[0])
    Pres  = float(f.readline().split(':')[0])
    u_inf = float(f.readline().split(':')[0])*1E+0
    Power = float(f.readline().split(':')[0])*1E+0
    T_fi = 0.5*(T_si+T_fl)
    print("T film = ", T_fi)
    
#    x = np.linspace(xmin,xmax,npt)

    N_Body = 0
    Cond = []
    Cp   = []
    Visc = []
    Pran = []
    Alfa = []
    Beta = []
    Body = []
    Air_OnOff   = int(f.readline().split(':')[0])
    if Air_OnOff == 1 :
        N_Body += 1
        MyData_Air = pd.read_excel('TermoPhys.xlsx',sheet_name=0).set_index('T')
        if T_fi > MyData_Air.index[len(MyData_Air)-1]:
            print("Temperature too high, choose one lower than :", MyData_Air.index[len(MyData_Air)-1])
            sys.exit(0) 
        Body.append('Air_T_fl')
        for i in range(len(MyData_Air)-1):
            if T_fi >= MyData_Air.index[i] and T_fi < MyData_Air.index[i+1]:
                a = (MyData_Air.iloc[i+1, 1] - MyData_Air.iloc[i, 1])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 1] - a*MyData_Air.index[i] 
                cp_value = (a*T_fi + b)*1.E+3
                Cp.append(cp_value)
                a = (MyData_Air.iloc[i+1, 3] - MyData_Air.iloc[i, 3])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 3] - a*MyData_Air.index[i] 
                visc_value = (a*T_fi + b)*1.E-6
                Visc.append(visc_value)
                a = (MyData_Air.iloc[i+1, 4] - MyData_Air.iloc[i, 4])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 4] - a*MyData_Air.index[i] 
                cond_value = (a*T_fi + b)*1.E-3
                Cond.append(cond_value)
                a = (MyData_Air.iloc[i+1, 5] - MyData_Air.iloc[i, 5])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 5] - a*MyData_Air.index[i] 
                alfa_value = (a*T_fi + b)*1.E-6
                Alfa.append(alfa_value)
                a = (MyData_Air.iloc[i+1, 6] - MyData_Air.iloc[i, 6])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 6] - a*MyData_Air.index[i] 
                pran_value = a*T_fi + b
                Pran.append(pran_value)
                a = (MyData_Air.iloc[i+1, 7] - MyData_Air.iloc[i, 7])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                b = MyData_Air.iloc[i, 7] - a*MyData_Air.index[i] 
                beta_value = (a*T_fi + b)*1.E-3
                Beta.append(beta_value)
                RayD = 9.81*beta_value*(T_si-T_fl)*Diam**3/visc_value/alfa_value
#                print(cp_value,visc_value,cond_value,pran_value,beta_value,alfa_value,RayD)
        
    Water_OnOff   = int(f.readline().split(':')[0])
    if Water_OnOff == 1 :
        N_Body  += 1
        MyData_Water = pd.read_excel('TermoPhys.xlsx',sheet_name=1).set_index('T')
        if T_fi > MyData_Water.index[len(MyData_Water)-1]:
            print("@@@@@@@@@@@@@@@@---------------&&&&&&&&&&&&&&&&&&&&&&")
            print("Temperature too high, choose one lower than :", MyData_Water.index[len(MyData_Water)-1])
            sys.exit(0) 
        Body.append('Water_T_fl')
        for i in range(len(MyData_Water)-1):
            if T_fi >= MyData_Water.index[i] and T_fi < MyData_Water.index[i+1]:
                a = (MyData_Water.iloc[i+1, 0] - MyData_Water.iloc[i, 0])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                b = MyData_Water.iloc[i, 0] - a*MyData_Water.index[i] 
                cp_value = (a*T_fi + b)*1.E+3
                Cp.append(cp_value)
                a = (MyData_Water.iloc[i+1, 2] - MyData_Water.iloc[i, 2])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                b = MyData_Water.iloc[i, 2] - a*MyData_Water.index[i] 
                visc_value = (a*T_fi + b)*1.E-6/997
                Visc.append(visc_value)
                a = (MyData_Water.iloc[i+1, 4] - MyData_Water.iloc[i, 4])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                b = MyData_Water.iloc[i, 4] - a*MyData_Water.index[i] 
                cond_value = (a*T_fi + b)*1.E-3
                Cond.append(cond_value)
                a = (MyData_Water.iloc[i+1, 6] - MyData_Water.iloc[i, 6])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                b = MyData_Water.iloc[i, 6] - a*MyData_Water.index[i] 
                pran_value = a*T_fi + b
                Pran.append(pran_value)
                a = (MyData_Water.iloc[i+1, 8] - MyData_Water.iloc[i, 8])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                b = MyData_Water.iloc[i, 8] - a*MyData_Water.index[i] 
                beta_value = (a*T_fi + b)*1.E-6
                Beta.append(beta_value)
                alfa = cond_value/cp_value/997
                RayD = 9.81*beta_value*(T_si-T_fl)*Diam**3/visc_value/alfa
#                print(cp_value,visc_value,cond_value,pran_value,beta_value,alfa,RayD)
        
    line=f.readline()   
    s = f.readline().split()
    if s[0][0] != '#' :
        s1,s2 = [float(s[0])*1.E-6, float(s[1])]
        N_Body  += 1
        Cond.append(0)
        Cp.append(0)
        Visc.append(s1)
        Pran.append(s2)
        Body.append(s[3])
    line=f.readline()
    while True:
        s = line.split()
        if s[0][0] != '#' :
            s1,s2 = [float(s[0])*1.E-6, float(s[1])]
            N_Body += 1
            Cond.append(0)
            Cp.append(0)
            Visc.append(s1)
            Pran.append(s2)
            Body.append(s[3])
        line=f.readline()
        if not line: break
    f.close()    

