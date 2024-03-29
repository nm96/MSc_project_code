# Lubrication Script
# ------------------

import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
import matplotlib.pyplot as plt
import time

from solutions import *
from models import *

t0 = time.time()
fn = 0  # Initialize figure number for plotting

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

# Define parameter values:

eps = 0.2 # Rotor eccentricity
Om = 4.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.5 # Damping coefficient
k = 10 # Stiffness coefficient
h = 1 # Gap width
mu = 1*10**-7 # Viscosity
b = 0.1 # Bearing length
R2 = 100 # Radius

Om_nat = (k/m)**0.5 # Shaft natural frequency

# Define the model:
model = (NHSommerfeld,(Om,h,mu,b,R2))
params = (eps,Om,m,c,k,model)

tspan = (0,2**10)    
N = tspan[1]*2**6
tt = np.linspace(*tspan,N)
X0 = [0.0001,0,0,0]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')

# Plot solution in stationary frame:
fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(*solxy(sol))
ax.plot(h*cos(np.linspace(0,2*pi,1000)),h*sin(np.linspace(0,2*pi,1000)),c='r') 
ax.set_aspect('equal')
ax.set_title("Solution trajectory in the stationary frame")

# Plot solution in corotating frame:
fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(*rotsolxy(sol,Om))
ax.plot(h*cos(np.linspace(0,2*pi,1000)),h*sin(np.linspace(0,2*pi,1000)),c='r') 
ax.set_aspect('equal')
ax.set_title("Solution trajectory in the co-rotating frame")


# Plot spectrum:
fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.axvline(Om_nat,ls='--',c='g')
ax.axvline(Om,ls='--',c='r')
ax.plot(*transformed(sol),c='k')
ax.set_title("Log-fft spectrum for a solution with the NHS lubrication model")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Log(fft(sol))")
ax.grid("on")


tf = time.time()
print("T = {:.2f}s".format(tf-t0))
plt.show()
