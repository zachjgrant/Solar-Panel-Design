#!/usr/bin/env python

#Script for testing function files

from func import n, delta, I_o, I_T, S
import numpy as np
import pandas as pd

##Function handles
# n(mon,day)
# delta(n)
# I_o(phi,delta,n)
# I_T(I_o,kt,rho_g,phi,delta,beta,gamma,n)
# S(I_o,kt,rho_g,phi,delta,beta,gamma,n,N,L,alpha_n)

mon = 7;
day = 14;
phi = 43.9; #deg
kt = .4;
beta = 30; #deg
gamma = 0;
rho_g = .4;
N = 2;
L = 0.003; #m
alpha_n = .93;

n = n(mon,day);
delta = delta(n)
I_o = I_o(phi,delta,n);
I_T = I_T(I_o,kt,rho_g,phi,delta,beta,gamma,n);
S = S(I_o,kt,rho_g,phi,delta,beta,gamma,n,N,L,alpha_n);

a = np.array([1,2,3,4,5,6]);
b = 1*5 - a*2;

print(I_T);
print(S);

