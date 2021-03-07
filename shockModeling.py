# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 15:14:55 2021

@author: Tagew
"""

import numpy as np
from scipy.integrate import odeint
import math
import matplotlib.pyplot as plt

# function that returns dz/dt

#Properties of argon
ma = 2
M = 39.948 #g/mol
cp= 0.5203
cv = 0.3122
gamma = cp/cv 
mu = 1 #viscosity
rho = 1 #Density
R = 8.314 #gas constant
T = 530
v1 = ma*math.sqrt(gamma*R*T/M)

def alfa(gamma, ma):
    f = (gamma-1)/(gamma+1) + 2/(ma**2*(gamma +1))
    return(f)
def beta(gamma):
    f = (9/8)*(gamma+1)*math.sqrt(math.pi/(8*gamma))
    return f
def lamb(T):
    f = (3*mu/rho)*math.sqrt(math.pi*M/8*R*T)
    return f


# function that returns dy/dt
def model(v,t):
    #dvdx = beta(gamma)*ma*(v-1)*(v-alfa(gamma,ma))/v
    # Forenklet
    dvdx = (v-1)*(v-2)/v
    return dvdx

# initial condition
y0 = 1

# time points
t = np.linspace(0,1000,10000)

# solve ODE
y = odeint(model,y0,t)

# plot results
plt.plot(t,y)
plt.xlabel('x*')
plt.ylabel('v*(x*)')
plt.show()