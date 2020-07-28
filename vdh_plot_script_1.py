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

# Define parameter values:
#eps = 0.14 # Rotor eccentricity
Om = 3.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.05 # Damping coefficient
k = 1 # Stiffness coefficient
h = 1 # Gap width
k_c = 50 # Stator stiffness parameter for VdH model

Om_nat = (k/m)**0.5 # Shaft natural frequency
Om_nat_c = (k_c/m)**0.5 # 'In-contact' natural frequency


# Define the model:
model = (VdH,(h,k_c))
    
tspan = (0,2**12)    
N = tspan[1]*2**6
tt = np.linspace(*tspan,N)
X0 = [0.01,0,0,0]

for eps in [0.12, 0.2, 0.30, 0.35, 0.38]:
    print(eps)
    params = (eps,Om,m,c,k,model)
    sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')
    # Plot spectrum:
    fn += 1; fig = plt.figure(fn)
    ax = fig.add_subplot(121)
    ax.axvline(Om_nat,ls='--',c='g')
    ax.axvline(Om,ls='--',c='r')
    ax.plot(*transformed(sol),c='k')
    ax.set_title(r"""$\varepsilon = {:.2f}$""".format(eps))
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Log(fft(sol))")
    ax.grid("on")
    ax = fig.add_subplot(122)
    ax.plot(*rotsolxy(sol,Om))
    ax.plot(h*cos(np.linspace(0,2*pi,1000)),h*sin(np.linspace(0,2*pi,1000)),c='r') 
    ax.set_aspect('equal')


tf = time.time()
print("T = {:.2f}s".format(tf-t0))

plt.show()
