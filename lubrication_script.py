import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
import matplotlib.pyplot as plt

from solutions import *
from models import *

fn = 0  # Initialize figure number for plotting

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

# Define parameter values:

eps = 0.03 # Rotor eccentricity
Om = 4.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.01 # Damping coefficient
k = 10 # Stiffness coefficient
h = 0.16 # Gap width
mu = 10**-5 # Viscosity
b = 1 # Bearing length
R2 = 1 # Radius

# Parameter values from Flores 2009

#h = 0.2*10**-3 # Clearance (m)
#b = 40*10**-3 # Bearing length (m)
#R2 = 9.8*10**-3 # Journal Radius (m)
#mu = 4 # Oil Viscosity (Pascals) - at 40C 
#m = 0.13 # Journal mass (kg) - big ?? on this one..
#k = 10 # ??
#c = 1 # ??
#eps = 0.3
#Om = 4.1
#

Om_nat = (k/m)**0.5 # Shaft natural frequency

# Define the model:
model = (NHSommerfeld,(Om,h,mu,b,R2))
params = (eps,Om,m,c,k,h,model)

tspan = (0,2**10)    
N = tspan[1]*2**6
tt = np.linspace(*tspan,N)
X0 = [0.0001,0,0,0]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')

# Plot solution in stationary frame:
fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(*solxy2(sol))
ax.plot(h*cos(np.linspace(0,2*pi,1000)),h*sin(np.linspace(0,2*pi,1000)),c='r') 
ax.set_aspect('equal')
ax.set_title("Solution trajectory in the stationary frame")

# Plot solution in corotating frame:
fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(*rotsolxy2(sol,Om))
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

plt.show()
