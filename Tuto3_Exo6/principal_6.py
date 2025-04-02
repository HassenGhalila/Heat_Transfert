#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 18:50:26 2023

@author: hassenghalila

COMMENTS:  
    
"""
import affect as af
import sortie as so
import operation as op
import pandas as pd


af.Affect()

nn = 0
if af.RayD < af.Ra :
    op.ChurchillChu()
    if af.Air_OnOff == 1 :
        while op.Ddt >= 1.E-6 :
    #        print(af.T_si)
            T_fi = 0.5*(op.T_s + af.T_fl)
            print("T film = ", T_fi)
            MyData_Air = pd.read_excel('TermoPhys.xlsx',sheet_name=0).set_index('T')
            for i in range(len(MyData_Air)-1):
                if T_fi >= MyData_Air.index[i] and T_fi < MyData_Air.index[i+1]:
                    a = (MyData_Air.iloc[i+1, 1] - MyData_Air.iloc[i, 1])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                    b = MyData_Air.iloc[i, 1] - a*MyData_Air.index[i] 
                    cp_value = (a*T_fi + b)*1.E+3
                    af.Cp[0] = cp_value
                    a = (MyData_Air.iloc[i+1, 3] - MyData_Air.iloc[i, 3])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                    b = MyData_Air.iloc[i, 3] - a*MyData_Air.index[i] 
                    visc_value = (a*T_fi + b)*1.E-6
                    af.Visc[0] = visc_value
                    a = (MyData_Air.iloc[i+1, 4] - MyData_Air.iloc[i, 4])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                    b = MyData_Air.iloc[i, 4] - a*MyData_Air.index[i] 
                    cond_value = (a*T_fi + b)*1.E-3
                    af.Cond[0] = cond_value
                    a = (MyData_Air.iloc[i+1, 5] - MyData_Air.iloc[i, 5])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                    b = MyData_Air.iloc[i, 5] - a*MyData_Air.index[i] 
                    alfa_value = (a*T_fi + b)*1.E-6
                    af.Alfa[0] = alfa_value
                    a = (MyData_Air.iloc[i+1, 6] - MyData_Air.iloc[i, 6])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                    b = MyData_Air.iloc[i, 6] - a*MyData_Air.index[i] 
                    pran_value = a*T_fi + b
                    af.Pran[0] = pran_value
                    a = (MyData_Air.iloc[i+1, 7] - MyData_Air.iloc[i, 7])/(MyData_Air.index[i+1]-MyData_Air.index[i])
                    b = MyData_Air.iloc[i, 7] - a*MyData_Air.index[i] 
                    beta_value = (a*T_fi + b)*1.E-3
                    af.Beta[0] = beta_value
                    RayD = 9.81*beta_value*(op.T_s-af.T_fl)*af.Diam**3/visc_value/alfa_value
#                    print(cp_value,visc_value,cond_value,pran_value,beta_value,alfa_value,RayD)
            nn += 1
            if nn > 20 : break
            op.ChurchillChu()
        
    if af.Water_OnOff == 1 :
        while op.Ddt >= 1.E-6 :
    #        print(af.T_si)
            T_fi = 0.5*(op.T_s + af.T_fl)
            print("T film = ", T_fi)
            MyData_Water = pd.read_excel('TermoPhys.xlsx',sheet_name=1).set_index('T')
            for i in range(len(MyData_Water)-1):
                if T_fi >= MyData_Water.index[i] and T_fi < MyData_Water.index[i+1]:
                    a = (MyData_Water.iloc[i+1, 0] - MyData_Water.iloc[i, 0])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                    b = MyData_Water.iloc[i, 0] - a*MyData_Water.index[i] 
                    cp_value = (a*T_fi + b)*1.E+3
                    af.Cp[0] = cp_value
                    a = (MyData_Water.iloc[i+1, 2] - MyData_Water.iloc[i, 2])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                    b = MyData_Water.iloc[i, 2] - a*MyData_Water.index[i] 
                    visc_value = (a*T_fi + b)*1.E-6/997
                    af.Visc[0] = visc_value
                    a = (MyData_Water.iloc[i+1, 4] - MyData_Water.iloc[i, 4])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                    b = MyData_Water.iloc[i, 4] - a*MyData_Water.index[i] 
                    cond_value = (a*T_fi + b)*1.E-3
                    af.Cond[0] = cond_value
                    a = (MyData_Water.iloc[i+1, 6] - MyData_Water.iloc[i, 6])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                    b = MyData_Water.iloc[i, 6] - a*MyData_Water.index[i] 
                    pran_value = a*T_fi + b
                    af.Pran[0] = pran_value
                    a = (MyData_Water.iloc[i+1, 8] - MyData_Water.iloc[i, 8])/(MyData_Water.index[i+1]-MyData_Water.index[i])
                    b = MyData_Water.iloc[i, 8] - a*MyData_Water.index[i] 
                    beta_value = (a*T_fi + b)*1.E-6
                    af.Beta[0]= beta_value
                    alfa = cond_value/cp_value/997
                    RayD = 9.81*beta_value*(op.T_s-af.T_fl)*af.Diam**3/visc_value/alfa
#                    print(cp_value,visc_value,cond_value,pran_value,beta_value,alfa,RayD)
            nn += 1
            if nn > 10 : break
            op.ChurchillChu()
else :
    print('Hi')
          
#so.ecrire()
