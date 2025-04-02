#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 18 21 18:50:26 2018

@author: hassenghalila

L'operation de trie est faite en utilisant les outils python
itemgetter permet de choisir l'element de reference pour le tri
"""
import affect as af
import numpy as np


def ChurchillChu():
    global Nu_D,h_bar,Ddt,T_s
    
    Nu_D  = np.zeros((af.N_Body))
    h_bar = np.zeros((af.N_Body))
    for i in range(af.N_Body):
        Nu_D[i] = (0.6 + 0.387*af.RayD**(1/6)/(1+(0.559/af.Pran[i])**(9/16))**(8/27))**2
        h_bar[i]  = af.Cond[i]*Nu_D[i]/af.Diam
        T_s = af.T_fl + af.Power/af.Diam/af.xmax/h_bar[i]/np.pi

        print("Ts = ", T_s)

    Ddt = abs(T_s - af.T_si)/af.T_fi
    af.T_si = T_s
    print("Erreur = ", Ddt)
    print("-------------------------- ")
