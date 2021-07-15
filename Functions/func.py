#!/usr/bin/env python

#Function file

#Input variables:
# mon = month of the year (1-12)
# day = day of the month
# phi = latitude
# I_o = extraterrestrial radiation on a horizontal surface for a given day of a given month at a given latitude
# kt = monthly average clearness index for location
# beta = angle of collector
# gamma = surface azimuth angle of collector
# rho_g = ground reflectance
# n = day of the year
# N = number of covers
# L = thickness of glass cover(s)
# alpha_n = absorptance of the plate at normal incidence

import numpy as np
import pandas as pd
from scipy.optimize import fsolve

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
def I_T(I_o,kt,rho_g,phi,delta,beta,gamma,n):

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

	#Initialize array to find where angle of incidence is 90 deg
	#o = np.linspace(0,180,1e6);

	#Initialize step array
	#s = np.linspace(0,1e6-1,1e6);
		
	#Calculate the sunset hour angle for the day
	if -np.tan(phi*dtr)*np.tan(delta*dtr) > 1 or -np.tan(phi*dtr)*np.tan(delta*dtr) < -1:
		omega = np.zeros(25);
		omega_s = 0;
	else:
		omega_s = np.degrees(np.arccos(-np.tan(phi*dtr)*np.tan(delta*dtr)));
		
	#Calculate when angle of incidence in 90 deg
	def f(x):
		return np.sin(delta*dtr)*np.sin(phi*dtr)*np.cos(beta*dtr) - np.sin(delta*dtr)*np.cos(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr) + np.cos(delta*dtr)*np.cos(phi*dtr)*np.cos(beta*dtr)*np.cos(x*dtr) + np.cos(delta*dtr)*np.sin(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr)*np.cos(x*dtr) + np.cos(delta*dtr)*np.sin(beta*dtr)*np.sin(gamma*dtr)*np.sin(x*dtr);
	starting_guess = omega_s;
	omega_0 = fsolve(f,starting_guess);

	#Correct the hour angle matrix to account for night in calculations of I_o
	omega = np.where(omega > omega_0, omega_0, omega);
	omega = np.where(omega < -omega_0, -omega_0, omega);

	#Initialize step array
	hour = np.linspace(0,23,24, dtype=int);

	#Create new hour angle matrix for calculations
	omega_n = np.zeros(24);
	for i in hour:
		omega_n[i] = (omega[i]+omega[i+1])/2;

	#Calculate R_b
	cos0 = np.zeros(24);
	cos0z = np.zeros(24);
	for i in hour:
		if omega_n[i] == omega_0 or omega_n[i] == -omega_0:
			cos0[i] = 0;
			cos0z[i] = 1e9;
		else:
			cos0[i] = np.sin(delta*dtr)*np.sin(phi*dtr)*np.cos(beta*dtr) - np.sin(delta*dtr)*np.cos(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr) + np.cos(delta*dtr)*np.cos(phi*dtr)*np.cos(beta*dtr)*np.cos(omega_n[i]*dtr) + np.cos(delta*dtr)*np.sin(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr)*np.cos(omega_n[i]*dtr) + np.cos(delta*dtr)*np.sin(beta*dtr)*np.sin(gamma*dtr)*np.sin(omega_n[i]*dtr);
			cos0z[i] = np.sin(delta*dtr)*np.sin(phi*dtr) + np.cos(delta*dtr)*np.cos(phi*dtr)*np.cos(omega_n[i]*dtr);
	R_b = cos0/cos0z;

	#Calculate I_T
	I_T = I_b*R_b + I_d*((1+np.cos(beta*dtr))/2) + I*rho_g*((1-np.cos(beta*dtr))/2);

	return I_T;

#Function for calculating the radiation absorbed by the collector
def S(I_o,kt,rho_g,phi,delta,beta,gamma,n,N,L,alpha_n):

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

	#Initialize array to find where angle of incidence is 90 deg
	#o = np.linspace(0,180,1e6);

	#Initialize step array
	#s = np.linspace(0,1e6-1,1e6);
		
	#Calculate the sunset hour angle for the day
	if -np.tan(phi*dtr)*np.tan(delta*dtr) > 1 or -np.tan(phi*dtr)*np.tan(delta*dtr) < -1:
		omega = np.zeros(25);
		omega_s = 0;
	else:
		omega_s = np.degrees(np.arccos(-np.tan(phi*dtr)*np.tan(delta*dtr)));
		
	#Calculate when angle of incidence in 90 deg
	def f(x):
		return np.sin(delta*dtr)*np.sin(phi*dtr)*np.cos(beta*dtr) - np.sin(delta*dtr)*np.cos(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr) + np.cos(delta*dtr)*np.cos(phi*dtr)*np.cos(beta*dtr)*np.cos(x*dtr) + np.cos(delta*dtr)*np.sin(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr)*np.cos(x*dtr) + np.cos(delta*dtr)*np.sin(beta*dtr)*np.sin(gamma*dtr)*np.sin(x*dtr);
	starting_guess = omega_s;
	omega_0 = fsolve(f,starting_guess);

	#Correct the hour angle matrix to account for night in calculations of I_o
	omega = np.where(omega > omega_0, omega_0, omega);
	omega = np.where(omega < -omega_0, -omega_0, omega);

	#Initialize step array
	hour = np.linspace(0,23,24, dtype=int);

	#Create new hour angle matrix for calculations
	omega_n = np.zeros(24);
	for i in hour:
		omega_n[i] = (omega[i]+omega[i+1])/2;

	#Calculate R_b
	cos0 = np.zeros(24);
	cos0z = np.zeros(24);
	for i in hour:
		if omega_n[i] == omega_0 or omega_n[i] == -omega_0:
			cos0[i] = 0;
			cos0z[i] = 1e9;
		else:
			cos0[i] = np.sin(delta*dtr)*np.sin(phi*dtr)*np.cos(beta*dtr) - np.sin(delta*dtr)*np.cos(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr) + np.cos(delta*dtr)*np.cos(phi*dtr)*np.cos(beta*dtr)*np.cos(omega_n[i]*dtr) + np.cos(delta*dtr)*np.sin(phi*dtr)*np.sin(beta*dtr)*np.cos(gamma*dtr)*np.cos(omega_n[i]*dtr) + np.cos(delta*dtr)*np.sin(beta*dtr)*np.sin(gamma*dtr)*np.sin(omega_n[i]*dtr);
			cos0z[i] = np.sin(delta*dtr)*np.sin(phi*dtr) + np.cos(delta*dtr)*np.cos(phi*dtr)*np.cos(omega_n[i]*dtr);
	R_b = cos0/cos0z;
	
	## BEAM RADIATION
	#Calculate angle of incidence after light moves through glass cover(s)
	theta_1_b = np.arccos(cos0); #radians
	n1_b = 1;
	n2_b = 1.526;
	theta_2_b = np.arcsin((n1_b/n2_b)*np.sin(theta_1_b)); #radians
	
	#Calculate reflectance values
	r_perp_b = ((np.sin(theta_2_b-theta_1_b))**2) / ((np.sin(theta_2_b+theta_1_b))**2);
	r_par_b = ((np.tan(theta_2_b-theta_1_b))**2) / ((np.tan(theta_2_b+theta_1_b))**2);
	
	#Calculate transmittance for polarization
	tau_r_b = .5 * (((1-r_par_b)/(1+(2*N-1)*r_par_b)) + ((1-r_perp_b)/(1+(2*N-1)*r_perp_b)));

	#Calculate the transmittance for absorption
	K = 4; #extinction coefficient for clear glass
	tau_a_b = np.exp((-K*L*N)/np.cos(theta_2_b));

	#Calculate total transmittance
	tau_b = tau_a_b * tau_r_b;

	#Calculate absorptance
	alpha_b = alpha_n * (1 - 1.5879e-3*theta_1_b + 2.7314e-4*theta_1_b**2 - 2.3026e-5*theta_1_b**3 + 9.0244e-7*theta_1_b**4 - 1.8000e-7*theta_1_b**5 + 1.7734e-10*theta_1_b**6 - 6.9937e-13*theta_1_b**7);

	#Calculate transmission-absorption constant
	tau_alpha_b = 1.01 * tau_b * alpha_b;

	## DIFFUSE RADIATION
	#Calculate angle of incidence after light moves through glass cover(s)
	theta_1_d = (59.7 - 0.1388*beta + 0.001497*beta**2) * dtr; #radians
	n1_d = 1;
	n2_d = 1.526;
	theta_2_d = np.arcsin((n1_d/n2_d)*np.sin(theta_1_d)); #radians
	
	#Calculate reflectance values
	r_perp_d = ((np.sin(theta_2_d-theta_1_d))**2) / ((np.sin(theta_2_d+theta_1_d))**2);
	r_par_d = ((np.tan(theta_2_d-theta_1_d))**2) / ((np.tan(theta_2_d+theta_1_d))**2);
	
	#Calculate transmittance for polarization
	tau_r_d = .5 * (((1-r_par_d)/(1+(2*N-1)*r_par_d)) + ((1-r_perp_d)/(1+(2*N-1)*r_perp_d)));

	#Calculate the transmittance for absorption
	tau_a_d = np.exp((-K*L*N)/np.cos(theta_2_d));

	#Calculate total transmittance
	tau_d = tau_a_d * tau_r_d;

	#Calculate absorptance
	alpha_d = alpha_n * (1 - 1.5879e-3*theta_1_d + 2.7314e-4*theta_1_d**2 - 2.3026e-5*theta_1_d**3 + 9.0244e-7*theta_1_d**4 - 1.8000e-7*theta_1_d**5 + 1.7734e-10*theta_1_d**6 - 6.9937e-13*theta_1_d**7);

	#Calculate transmission-absorption constant
	tau_alpha_d = 1.01 * tau_d * alpha_d;

	## GROUND REFLECTED RADIATION
	#Calculate angle of incidence after light moves through glass cover(s)
	theta_1_g = (90 - 0.5788*beta + 0.002693*beta**2) * dtr; #radians
	n1_g = 1;
	n2_g = 1.526;
	theta_2_g = np.arcsin((n1_d/n2_g)*np.sin(theta_1_g)); #radians
	
	#Calculate reflectance values
	r_perp_g = ((np.sin(theta_2_g-theta_1_g))**2) / ((np.sin(theta_2_g+theta_1_g))**2);
	r_par_g = ((np.tan(theta_2_g-theta_1_g))**2) / ((np.tan(theta_2_g+theta_1_g))**2);
	
	#Calculate transmittance for polarization
	tau_r_g = .5 * (((1-r_par_g)/(1+(2*N-1)*r_par_g)) + ((1-r_perp_g)/(1+(2*N-1)*r_perp_d)));

	#Calculate the transmittance for absorption
	tau_a_g = np.exp((-K*L*N)/np.cos(theta_2_g));

	#Calculate total transmittance
	tau_g = tau_a_g * tau_r_g;

	#Calculate absorptance
	alpha_g = alpha_n * (1 - 1.5879e-3*theta_1_g + 2.7314e-4*theta_1_g**2 - 2.3026e-5*theta_1_g**3 + 9.0244e-7*theta_1_g**4 - 1.8000e-7*theta_1_g**5 + 1.7734e-10*theta_1_g**6 - 6.9937e-13*theta_1_g**7);

	#Calculate transmission-absorption constant
	tau_alpha_g = 1.01 * tau_g * alpha_g;

	#Calculate S
	S = I_b*R_b*tau_alpha_b + I_d*tau_alpha_d*((1+np.cos(beta*dtr))/2) + I*rho_g*tau_alpha_g*((1-np.cos(beta*dtr))/2);

	return S;




