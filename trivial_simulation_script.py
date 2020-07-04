import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
import matplotlib.pyplot as plt

from solutions import *
from models import *

fn = 0  # Initialize figure number for plotting

# Define parameter values:
eps = 0.0525 # Rotor eccentricity
Om = 6.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.3 # Damping coefficient
k = 10 # Stiffness coefficient
h = 0.05 # Gap width

Om_nat = (k/m)**0.5 # Shaft natural frequency


# Trivial model
# -------------

model = (simple,(0,)) 

params = (eps,Om,m,c,k,h,model) # Package parameters into a tuple
    
# Integrate the ODE system over a given time span with given initial
# conditions:

tspan = (0,2**12)    
N = tspan[1]*2**6
tt = np.linspace(*tspan,N)
X0 = [0.1,0,-0.2,0]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params)

fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(sol.t,sol.y[2])

fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(*solxy(sol))
ax.set_title("Solution trajectory in the stationary frame")

fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(*rotsolxy(sol,Om))
ax.set_title("Solution trajectory in the rotating frame")

fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])

ax.axvline(Om_nat,ls='--',c='g')
ax.axvline(Om,ls='--',c='r')
ax.plot(*transformed(sol))
ax.set_title("Log-fft spectrum for a solution with the trivial model")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Log(fft(sol))")
ax.grid("on")

plt.show()
