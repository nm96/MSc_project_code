import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
import matplotlib.pyplot as plt
import time
from solutions import *
from models import *
from matplotlib import rc
import sys

# Set a nice font for the plots when on Linux
if sys.platform == 'linux':
    rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    rc('text', usetex=True)

t0 = time.time()
fn = 0  # Initialize figure number for plotting

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

# Integration limits
tspan = (0,2**12)    
N = tspan[1]*2**6
tt = np.linspace(*tspan,N)
X0 = [0.0001,0,0,0]

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
Om = 0.3
c = 0.10

fn += 1; fig = plt.figure(fn,figsize=[6,6])
spn = 110
#for B in np.linspace(1.49,1.5,4):
for B in [1.5]:
    Om_nat = (k/m)**0.5 # Shaft natural frequency
    model = (NHSommerfeld2,(Om,h,B))
    params = (eps,Om,m,c,k,h,model)
    sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')
    # Plot spectrum:
    spn += 1
    ax = fig.add_subplot(spn)
    ax.axvline(Om_nat,ls='--',c='g')
    ax.axvline(Om,ls='--',c='r')
    ax.plot(*transformed(sol),c='k')
    ax.set_title(r"""$\beta$ = {:.4f}""".format(B))
    ax.set_ylabel("$\log|\mathcal{F}[X]|$")
    ax.set_xlabel("$\omega \ (s^{-1})$")
    ax.grid("on")

plt.tight_layout()

tf = time.time()
print("T = {:.2f}s".format(tf-t0))
plt.show()
