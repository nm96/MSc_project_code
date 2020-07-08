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

#Define parameter values:
#eps = 0.1 # Rotor eccentricity
#Om = 1.6 # Driving frequency
#m = 1 # Mass (per unit length)
#c = 0.01 # Damping coefficient
#k = 0 # Stiffness coefficient
#h = 0.2 # Gap width

# Define parameter values:
eps = 0.0525 # Rotor eccentricity
Om = 4.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.05 # Damping coefficient
k = 10 # Stiffness coefficient
h = 0.1 # Gap width

k_c = 50 # Stator stiffness parameter for models
alpha = 500 # for the Hua model
mu = 0.1 # Coefficient of friction

Om_nat = (k/m)**0.5 # Shaft natural frequency
Om_nat_c = (k_c/m)**0.5 # 'In-contact' natural frequency


# Trivial model
# -------------

# Define the model:
model = (simple,(0,)) 
params = (eps,Om,m,c,k,h,model) # Package parameters into a tuple
    
# Set conditions and time span, integrate

tspan = (0,2**10)    
N = tspan[1]*2**6
tt = np.linspace(*tspan,N)
X0 = [0.01,0,0,0]

sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')

# Plot spectrum:

fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(*solxy(sol))
ax.plot(h*cos(np.linspace(0,2*pi,1000)),h*sin(np.linspace(0,2*pi,1000)),c='r') 
ax.set_aspect('equal')
ax.set_title("Solution trajectory in the stationary frame")

fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])
ax.plot(*rotsolxy(sol,Om))
ax.plot(h*cos(np.linspace(0,2*pi,1000)),h*sin(np.linspace(0,2*pi,1000)),c='r') 
ax.set_aspect('equal')
ax.set_title("Solution trajectory in the co-rotating frame")

fn += 1; fig = plt.figure(fn); ax = fig.add_axes([.1,.1,.8,.8])

ax.axvline(Om_nat,ls='--',c='g')
ax.axvline(Om,ls='--',c='r')
ax.plot(*transformed(sol))
ax.set_title("Log-fft spectrum for a solution with the trivial model")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Log(fft(sol))")
ax.grid("on")


# Van der Heijden Model
# ---------------------

# Define the model:
model = (VdH,(h,k_c))
params = (eps,Om,m,c,k,h,model)
    
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
ax.axvline(Om_nat_c,ls='--',c='b')
ax.axvline(Om,ls='--',c='r')
ax.plot(*transformed(sol),c='k')
ax.set_title("Log-fft spectrum for a solution with the Van der Heijden model")
ax.set_xlabel("Frequency (Hz)")
ax.set_ylabel("Log(fft(sol))")
ax.grid("on")

plt.show()
