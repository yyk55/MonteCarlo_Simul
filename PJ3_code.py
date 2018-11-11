# -*- coding: utf-8 -*-
"""
Created on Sat Apr 28 21:17:05 2018

@author: na88555
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from datetime import date

# parameters
#############################################################
SS0 = 2614.45
RS0 = 1512.155
T = (date(2023, 9, 11) - date(2018, 4, 4)).days / 365
M = np.busday_count(date(2023, 6, 6), date(2023, 9,6))
dt = T / M
Pdate = date(2018, 4, 4)
Idate = date(2018, 4, 9)
Ssig = 0.20702
Rsig = 0.23259
r = 0.0274275
ro = 0.0274252
Sdiv = 0.02042
Rdiv = 0.011228
corr = 0.849
principal = 1000
num_simulation = 10000
#############################################################


def monte_carlo(SS0, RS0, T, M, dt, Pdate, Idate, Ssig, Rsig, r, ro, Sdiv, Rdiv, corr, principal, num_simulation):
	SS = np.zeros(M) 
	RS = np.zeros(M)
	SS[0] = SS0 
	RS[0] = RS0
	
	# Start Monte-carlo simulation
	avgs = np.zeros(num_simulation)
	avgr = np.zeros(num_simulation)
	VT = np.zeros(num_simulation)
	
	for num in range(num_simulation):
		for i in range(M-1):  # index paths 
			phi1 = np.random.standard_normal()
			phi2 = np.random.standard_normal()
			SS[i+1] = SS[i]*np.exp((r - Sdiv - 0.5 * Ssig**2)*dt + Ssig * np.sqrt(dt) * phi1)
			RS[i+1] = RS[i]*np.exp((r - Rdiv - 0.5 * Rsig**2)*dt + Rsig * corr * np.sqrt(dt) 
			* phi1 + Rsig * np.sqrt(1-corr**2) * np.sqrt(dt) * phi2)
		avgs[num] = np.average(SS) # final average index values
		avgr[num] = np.average(RS)
		pfs = avgs[num] / SS0
		pfr = avgr[num] / RS0
		if pfs > pfr:
			wpf = pfr
		else:
			wpf = pfs
		# payoffs
		if avgs[num] >= SS0 * 1.21 and avgr[num] >= RS0 * 1.21:
			VT[num] = principal + principal * ((wpf-1.21) * 3.34) + 415
		elif SS0 <= avgs[num] < SS0 * 1.21 or RS0 <= avgr[num] < RS0 * 1.21:
			VT[num] = principal + principal * ((wpf-1) * 1.50) + 100
		elif SS0*0.95 <= avgs[num] < SS0 or RS0*0.95 <= avgr[num] < RS0:
			VT[num] = principal + principal * ((wpf-0.95) * 2)
		elif avgs[num] < SS0*0.95 or avgr[num] < RS0*0.95:
			VT[num] = principal * wpf + 50
	
	# Final value of the product (with adjustment btwn Issue date and Pricing date)
	V0 = np.exp(ro * ((Idate-Pdate).days / 365)) * np.exp(-r * T) * 1/num_simulation * np.sum(VT)
	
	return V0

Final_V = monte_carlo(SS0, RS0, T, M, dt, Pdate, Idate, Ssig, Rsig, r, ro, Sdiv, Rdiv, corr, principal, num_simulation)

print(Final_V)		


Vs = np.zeros(50)
for i in range(50):
	fv = monte_carlo(SS0, RS0, T, M, dt, Pdate, Idate, Ssig, Rsig, r, ro, Sdiv, Rdiv, corr, principal, num_simulation)
	Vs[i] = fv		

fig, ax1 = plt.subplots()
ax1.plot(Vs, 'b')
ax1.set_xlabel('The number of trials with each 10,000 simulations')
ax1.set_ylabel('value of the product')

fig.tight_layout()
plt.show()	
	 
	