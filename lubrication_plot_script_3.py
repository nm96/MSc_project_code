import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
import matplotlib.pyplot as plt
import matplotlib.ticker
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
Om = 3*2**0.5

fn += 1; fig = plt.figure(fn,figsize=[12,7])
spn = 110
#for B in np.linspace(1.49,1.5,4):
for B in [0.5]:
    Om_nat = (k/m)**0.5 # Shaft natural frequency
    model = (NHSommerfeld2,(Om,h,B))
    params = (eps,Om,m,c,k,model)
    sol = solve_ivp(dXdt,tspan,X0,t_eval=tt,args=params,method='Radau')
    # Plot spectrum, label:
    spn += 1
    ax = fig.add_subplot(spn)
    ax.set_xlim([0,10])
    ax.axvline(Om_nat,ls='--',c='g',label=r"$\omega_{nat}$")
    ax.axvline(Om,ls='--',c='r',label=r"$\Omega$")
    om, P = PSD_nw(sol)
    ax.semilogy(om,P,c='gray',label=r"Spectrum without windowing")
    om, P = PSD(sol)
    ax.semilogy(om,P,c='k',label=r"Spectrum with Hanning window applied")
    locmaj = matplotlib.ticker.LogLocator(base=100,numticks=20) 
    ax.yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=100,subs=(0.2,0.4,0.6,0.8),numticks=50)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.grid()
    ax.set_title(r"""Effect of Windowing ($\beta = {:.2f}$)""".format(B))
    ax.set_ylabel("$P(\omega)$",rotation=0)
    ax.yaxis.labelpad = 20
    ax.set_xlabel("$\omega \ (s^{-1})$")
    ax.legend()

plt.tight_layout()

#fig.savefig("../plots/windowing_effects.eps")

tf = time.time()
print("T = {:.2f}s".format(tf-t0))
plt.show()
