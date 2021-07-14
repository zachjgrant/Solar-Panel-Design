#!/usr/bin/env python

#Function file

#Input variables:
# mon = month of the year (1-12)
# day = day of the month
# phi = latitude
# I_o = Extraterrestrial radiation on a horizontal surface for a given day of a given month at a given latitude
# kt = monthly average clearness index for location
# beta = angle of collector
# gamma = surface azimuth angle of collector

import numpy as np
import pandas as pd

#Helpful constants
pi = np.pi;
dtr = pi/180;
rtd = 1/dtr;

#Function for determining the day of the year
def n(mon,day):
	#Initialize array for days
	if mon == 1:
		#day = np.linspace(1,31,31); Will use for monthly version
		n = day;
		#l = 31; Will use for monthly version
	elif mon == 2:
		#day = np.linspace(32,59,28);
		n = 31 + day;
		#l = 28;
	elif mon == 3:
		#day = np.linspace(60,90,31);
		n = 59 + day;
		#l = 31;
	elif mon == 4:
		#day = np.linspace(91,120,30);
		n = 90 + day;
		#l = 30;
	elif mon == 5:
		#day = np.linspace(121,151,31);
		n = 120 + day;
		#l = 31;
	elif mon == 6:
		#day = np.linspace(152,181,30);
		n = 151 + day;
		#l = 30;
	elif mon == 7:
		#day = np.linspace(182,212,31);
		n = 181 + day;
		#l = 31;
	elif mon == 8:
		#day = np.linspace(213,243,31);
		n = 212 + day;
		#l = 31;
	elif mon == 9:
		#day = np.linspace(244,273,30);
		n = 243 + day;
		#l = 30;
	elif mon == 10:
		#day = np.linspace(274,304,31);
		n = 273 + day;
		#l = 31;
	elif mon == 11:
		#day = np.linspace(305,334,30);
		n = 304 + day;
		#l = 30;
	else:
		#day = np.linspace(335,365,31);
		n = 334 + day;
		#l = 31;

	return n;

#Function for calculating the declination
def delta(n):

	#Calculate declination for day
	delta = 23.45 * np.sin(2*pi*(284+n)/365);

	return delta;

#Function for calculating the extraterrestrial hourly radiation for a given day of a given month at a given latitude
def I_o(phi,delta,n):
	
	#Initializing array for hour angle
	omega = np.linspace(-180,180,25);
	
	#Initialize step array
	hour = np.linspace(0,23,24, dtype=int);
	
	#Initialize the storage array
	I_o = np.zeros(24);

	#Calculate the sunset hour angle for the day
	if -np.tan(phi*dtr)*np.tan(delta*dtr) > 1 or -np.tan(phi*dtr)*np.tan(delta*dtr) < -1:
		omega = np.zeros(25);
		omega_s = 0;
	else:
		omega_s = np.degrees(np.arccos(-np.tan(phi*dtr)*np.tan(delta*dtr)));
	
	#Correct the hour angle matrix to account for night in calculations of I_o
	omega = np.where(omega > omega_s, omega_s, omega);
	omega = np.where(omega < -omega_s, -omega_s, omega);	

	#Calculate I_o for each hour
	for i in hour:
		I_o[i] = (12*3600*1367/pi) * (1+0.033*np.cos(2*pi*n/365)) * ((np.cos(phi*dtr)*np.cos(delta*dtr)*np.sin((omega[i+1]-omega[i])*dtr))+((omega[i+1]-omega[i])*dtr*np.sin(phi*dtr)*np.sin(delta*dtr)));

	return I_o;

#Function for calculating the radiation on a tilted surface inside the earth's atmosphere
def I_T(I_o,kt,phi,delta,beta,gamma,n):

	#Calculate radiation on a horizontal surface inside the earth's atmosphere	
	I = kt * I_o;
	
	#Calculate diffuse radiation
	if kt <= 0.22:
		I_d = I * (1.0-0.09*kt);
	elif kt > 0.22 and kt <= 0.80:
		I_d = I * (0.9511-0.1604*kt+4.388*kt**2-16.638*kt**3+12.336*kt**4);
	else:
		I_d = I * 0.165;

	#Calculate beam radiation
	I_b = I - I_d;

	#Initialize array for hour angle
	omega = np.linspace(-180,180,25);

	#Calculate the sunset hour angle for the day
	if -np.tan(phi*dtr)*np.tan(delta*dtr) > 1 or -np.tan(phi*dtr)*np.tan(delta*dtr) < -1:
		omega = np.zeros(25);
		omega_s = 0;
	else:
		omega_s = np.degrees(np.arccos(-np.tan(phi*dtr)*np.tan(delta*dtr)));
	
	#Correct the hour angle matrix to account for night in calculations of I_o
	omega = np.where(omega > omega_s, omega_s, omega);
	omega = np.where(omega < -omega_s, -omega_s, omega);

	#Initialize step array
	hour = np.linspace(0,23,24, dtype=int);

	#Create new hour angle matrix for calculations
	omega_n = np.zeros(24);
	for i in hour:
		omega_n[i] = (omega[i]+omega[i+1])/2;

	#Calculate R_B
	cos0 = np.zeros(24);
	cos0z = np.zeros(24);
	for i in hour:
		if omega_n[i] == omega_s or omega_n[i] == -omega_s:
			cos0[i] = 0;
			cos0z[i] = 1e9;
		else:
			cos0[i] = np.sin(delta*dtr)*np.sin(phi*dtr)*np.cos(beta*dtr) - np.sin(delta*dtr)*np.cos(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr) + np.cos(delta*dtr)*np.cos(phi*dtr)*np.cos(beta*dtr)*np.cos(omega_n[i]*dtr) + np.cos(delta*dtr)*np.sin(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr)*np.cos(omega_n[i]*dtr) + np.cos(delta*dtr)*np.sin(beta*dtr)*np.sin(gamma*dtr)*np.sin(omega_n[i]*dtr);
			cos0z[i] = np.sin(delta*dtr)*np.sin(phi*dtr) + np.cos(delta*dtr)*np.cos(phi*dtr)*np.cos(omega_n[i]*dtr);
	R_B = cos0/cos0z;

	return R_B;






