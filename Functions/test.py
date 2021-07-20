#!/usr/bin/env python

#Script for testing function files

from func import n, delta, I_o, I_T, S, Q_u, eta
import numpy as np
import pandas as pd
from prettytable import PrettyTable

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
day = 19;
phi = 43.9; #deg
kt = 0.6;
beta = 15; #deg
gamma = 10;
rho_g = 0.4;
N = 2;
L = 0.003; #m
alpha_n = 0.93;
T_a = (np.array([61,61,61,61,61,61,63,65,67,70,73,75,76,77,78,79,75,70,69,68,67,65,65,62])-32)*5/9 + 273.15;
T_dp = (np.array([61,61,61,61,61,61,62,62,63,64,65,65,66,66,66,66,64,65,64,63,62,62,62,62])-32)*5/9 + 273.15;
T_fi = (120-32)*5/9+273.15; #120 deg F
L_p = 0.015; #m
e = np.array([0.95,0.88,0.88]);
W = 0.15; #m
D = 0.01; #m
D_i = 0.008; #m
delt = 0.0005; #m
k_p = 385; #W/mK
m_dot = 1.0; #kg/min
A_c = 2; #m^2

#Calculations
n = n(mon,day);
delta = delta(n)
I_o = I_o(phi,delta,n);
I_T = I_T(I_o,kt,rho_g,phi,delta,beta,gamma,n);
S = S(I_o,kt,rho_g,phi,delta,beta,gamma,n,N,L,alpha_n);
Q_u = Q_u(S,T_a,T_dp,T_fi,N,L_p,e,beta,W,D,D_i,delt,k_p,m_dot,A_c);
eta = eta(Q_u,I_T,A_c)*100;

#Create characteristics table and print
ch = PrettyTable(["Variable","Value"]);
ch.title = "Solar Collector Characteristics";
ch.add_row(["Date",str(mon)+"/"+str(day)]);
ch.add_row(["Latitude (\u03d5)",str(phi)+"\u00b0"]);
ch.add_row(["Collector Slope (\u03b2)",str(beta)+"\u00b0"]);
ch.add_row(["Collector Azimuth Angle (\u03b3)",str(gamma)+"\u00b0"]);
ch.add_row(["Daily Clearness Index (k\u209c)",str(kt)]);
ch.add_row(["Ground Reflectance (\u03c1)",str(rho_g)]);
ch.add_row(["# of Covers (N)", str(N)]);
ch.add_row(["Thickness of Cover (L)", str(L*1000) + " mm"]);
ch.add_row(["Absorptance of Plate at Normal Incidence (\u03b1\u2099)", str(alpha_n)]);
ch.add_row(["Fluid Inlet Temperature (T)",str(np.int32((T_fi-273.15)*9/5+32))+"\u00b0F"]);
ch.add_row(["Distance Between Plates (L\u209a)",str(L_p*1000)+" mm"]);
ch.add_row(["Emissivity of Plate (\u03b5\u209a)",str(e[0])]);
ch.add_row(["Distance Between Tubes (W)",str(W*1000)+" mm"]);
ch.add_row(["Outer Diamter of Tubes (D\u2092)",str(D*1000)+" mm"]);
ch.add_row(["Inner Diamter of Tubes (D)",str(D_i*1000)+" mm"]);
ch.add_row(["Thickness of Plate (\u03b4)",str(delt*1000)+" mm"]);
ch.add_row(["Thermal Conductivity of Plate (k\u209a)",str(k_p)+" W/mK"]);
ch.add_row(["Fluid Mass Flow Rate (m\u0307)",str(m_dot)+" kg/min"]);
ch.add_row(["Collector Area (A)",str(A_c)+" m\u00b2"]);


print(ch);

#Create output table and print
T_a_t = np.round((T_a-273.15)*(9/5)+32,1);
T_a_t = T_a_t.astype(int);
T_dp_t = np.round((T_dp-273.15)*(9/5)+32,1);
T_dp_t = T_dp_t.astype(int);
I_o_t = np.round(I_o/3600,1);
I_T_t = np.round(I_T/1e3,1);
S_t = np.round(S/1e3,1);
Q_u_t = np.round(Q_u/1e3,1);
eta_t = np.round(eta,1);
out = PrettyTable(["Hour","Ambient Temperature (F)", "Dew Point Temperature (F)","I_o (kJ)","I_T (kJ)","S (kJ)","Q_u (kJ)","Efficiency (%)"]);
out.title = "Energy Analysis";
out.padding_width = 0;
hour = ["0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12","12-13","13-14","14-15","15-16","16-17","17-18","18-19","19-20","20-21","21-22","22-23","23-24"];
for i in range(24):
	out.add_row([hour[i],T_a_t[i],T_dp_t[i],I_o_t[i],I_T_t[i],S_t[i],Q_u_t[i],eta_t[i]]);

print(out);

