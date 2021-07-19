#!/usr/bin/env python

#Script for testing function files

from func import n, delta, I_o, I_T, S, Q_u, eta
import numpy as np
import pandas as pd

##Function handles
# n(mon,day)
# delta(n)
# I_o(phi,delta,n)
# I_T(I_o,kt,rho_g,phi,delta,beta,gamma,n)
# S(I_o,kt,rho_g,phi,delta,beta,gamma,n,N,L,alpha_n)
# Q_u(S,T_a,T_dp,T_fi,N,L_p,e,beta,W,D,D_i,delt,k_p,m_dot,A_c)
# eta(Q_u,I_T,A_C)

##Variables
# mon = month of the year (1-12)
# day = day of the month
# phi = latitude (deg)
# I_o = extraterrestrial radiation on a horizontal surface  (J/m^2 per hour)
# kt = monthly average clearness index for location
# beta = angle of collector (deg)
# gamma = surface azimuth angle of collector (deg)
# rho_g = ground reflectance
# n = day of the year
# I_T = rasiation on the collector (J/m^2 per hour)
# N = number of covers
# L = thickness of glass cover(s) (m)
# alpha_n = absorptance of the plate at normal incidence
# S = radiation absorbed by the collector for the hour (J/m^2 per hour)
# T_a = ambient temperature as an array of 24 values (K)
# T_dp = dew point temperature as an array of 24 values (K)
# T_fi = fluid inlet temperature (K)
# L_p = distance between covers (m)
# e = emissivity of surfaces from plate to outside cover
# W = distance between tubes (m)
# D = outside diamter of tubes (m)
# D_i = inside diameter of tubes (m)
# delt = plate thickness (m)
# k_p = thermal conductivity of plate (W/mk) (Copper = 385 W/mK)
# m_dot = mass flow rate of fluid (kg/min)
# A_c = collector area (m^2)
# Q_u = useful energy gain from collector (J/hr)
# eta = collector efficiency

#Input
mon = 7;
day = 14;
phi = 43.9; #deg
kt = .6;
beta = 15; #deg
gamma = 10;
rho_g = .4;
N = 2;
L = 0.003; #m
alpha_n = .93;
T_a = np.array([20,20,19,19,18,18,17,18,20,21,23,24,25,29,29,29,29,30,28,27,25,23,22,21])+273;
T_dp = np.array([18,17,17,16,16,16,16,17,18,18,18,18,18,22,22,22,21,21,21,21,20,19,19,18])+273;
T_fi = 323; #50 deg C
L_p = 0.015;
e = np.array([.95,.88,.88]);
W = .15;
D = 0.01;
D_i = 0.008;
delt = 0.0005;
k_p = 385;
m_dot = 1;
A_c = 2;

#Calculations
n = n(mon,day);
delta = delta(n)
I_o = I_o(phi,delta,n);
I_T = I_T(I_o,kt,rho_g,phi,delta,beta,gamma,n);
S = S(I_o,kt,rho_g,phi,delta,beta,gamma,n,N,L,alpha_n);
Q_u = Q_u(S,T_a,T_dp,T_fi,N,L_p,e,beta,W,D,D_i,delt,k_p,m_dot,A_c);
eta = eta(Q_u,I_T,A_c);

print(I_o);
print(I_T);
print(S);
print(Q_u);
print(eta);

