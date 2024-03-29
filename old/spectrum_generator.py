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
tspan = (0,2**16)  
N = tspan[1]*2**4
tt = np.linspace(*tspan,N)
X0 = [0.01,0,0,0]

# Parameter values:
eps = 0.2 # Rotor eccentricity
Om = 4.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.5 # Damping coefficient
k = 10 # Stiffness coefficient
h = 1 # Gap width

B = 0.7
c = 0.4
Om = 0.3

# Solve and compute spectrum

model = (NHSommerfeld2,(Om,h,B))
params = (eps,Om,m,c,k,model)
sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')
om, P = PSD(sol)
filename = "beta_{:.2f}_Om_{:.2f}_spectrum.csv".format(B,Om)
np.savetxt(filename,np.array([om,P]).T)
print("Created", filename)

tf = time.time()
print("tf = {:.2f}s".format(tf-t0))
