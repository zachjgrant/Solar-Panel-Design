#!/usr/bin/env python

#Script for testing function files

from func import n, delta, I_o, I_T
import numpy as np
import pandas as pd

mon = 7;
day = 14;
phi = 43.9;
kt = .5;
beta = 30;
gamma = 0;

n = n(mon,day);
delta = delta(n)
I_o = I_o(phi,delta,n);
cos0 = I_T(I_o,kt,phi,delta,beta,gamma,n);

a = np.array([-2.8,-2.8,-2.8,-2,-1,0,1,2,2.8,2.8,2.8]);
b = np.array([0,1,2,3,4,5,6,7,8,9]);
c = np.zeros(10);
for i in b:
	c[i] = (a[i]+a[i+1])/2;

d = np.zeros(10);
k = 0;
for j in c:
	if j == -2.8 or j == 2.8:
		d[k] = 0;
	else:
		d[k] = 1;
	k = k+1;

print(cos0)

