# A script for integrating the Jeffcott equations using scipy's solve_ivp.

import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
import matplotlib.pyplot as plt

from solutions import *
from models import *

fn = 0  # Initialize figure number for plotting

# Define parameter values
eps = 0.0525 # Rotor eccentricity
Om = 3.5 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.05 # Damping coefficient
k = 10 # Stiffness coefficient
h = 0.03 # Gap width

model = (simple,(0,)) # Model-tuple for the trivial no-force model

params = (eps,Om,m,c,k,h,model) # Package parameters into a tuple
    
# Integrate the ODE system over a given time span with given initial conditions

N = 2**15
tspan = (0,1000)    
tt = np.linspace(*tspan,N)
X0 = [0.1,0,0,0]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params)

# Plot fft

max_freq = 10
freq_scale = np.arange(0,max_freq,4*pi/sol.t[-1])
print(2*pi/sol.t[1])

fn += 1; fig = plt.figure(fn)
ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(freq_scale,np.log(np.abs(rfft(sol.y[0][N//2:])))[:len(freq_scale)])
plt.grid("on")


# Integration with VdH model

# Define the model
k_c = 50 # Required model-specific parameter
model = (VdH,(h,k_c)) # This is the standard form for a model;
# A tuple containing a model function with the form f(X,mparams) and a tuple
# mparams of model parameters.

params = (eps,Om,m,c,k,h,model)
    
# Integrate the ODE system over a given time span with given initial conditions

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params)

# Plot fft

max_freq = 10
freq_scale = np.arange(0,max_freq,4*pi/sol.t[-1])
print(2*pi/sol.t[1])

fn += 1; fig = plt.figure(fn)
ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(freq_scale,np.log(np.abs(rfft(sol.y[0][N//2:])))[:len(freq_scale)])
plt.grid("on")

plt.show()
