import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
from solutions import *
from models import *
import time

t0 = time.time()
fn = 0  # Initialize figure number for plotting

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

# Integration limits
tspan = (0,2**14)  
N = tspan[1]*2**6
tt = np.linspace(*tspan,N)
X0 = [0.01,0,0,0]

# Standard parameter values:
eps = 0.2 # Rotor eccentricity
Om = 4.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.5 # Damping coefficient
k = 10 # Stiffness coefficient
h = 1 # Gap width
mu = 1*10**-7 # Viscosity
b = 0.1 # Bearing length
R2 = 100 # Radius

# Special values:
#c = 0.11
Om = 4.1

for B in [1.65]:
    model = (NHSommerfeld2,(Om,h,B))
    params = (eps,Om,m,c,k,model)
    sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')
    om, P = PSD(sol)
    np.savetxt("beta_{:.2f}_spectrum.csv".format(B),np.array([om,P]).T)

tf = time.time()
print("Time taken = {:.2f}s".format(tf-t0))
