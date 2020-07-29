import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib import rc
import sys

# Set a nice font for the plots when on Linux
if sys.platform == 'linux':
    rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    rc('text', usetex=True)

fn = 0  # Initialize figure number for plotting

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

fn += 1; fig = plt.figure(fn,figsize=[12,6])
spn = 320
#for B in np.linspace(1.49,1.5,4):
for B in [0.15, 0.80, 1.30, 1.35, 1.40, 1.65]:
    Om_nat = (k/m)**0.5 # Shaft natural frequency
    # Plot spectrum, label:
    spn += 1
    ax = fig.add_subplot(spn)
    ax.set_xlim([0,10])
    ax.axvline(Om_nat,ls='--',c='g',label=r"$\omega_{nat}$")
    ax.axvline(Om,ls='--',c='r',label=r"$\Omega$")
    om, P = np.loadtxt("beta_{:.2f}_spectrum.csv".format(B)).T
    ax.semilogy(om,P,c='k')
    locmaj = matplotlib.ticker.LogLocator(base=100,numticks=30) 
    ax.yaxis.set_major_locator(locmaj)
    locmin = matplotlib.ticker.LogLocator(base=100,subs=(0.2,0.4,0.6,0.8),numticks=50)
    ax.yaxis.set_minor_locator(locmin)
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.grid()
    ax.set_title(r"""$\beta = {:.2f}$""".format(B))
    ax.set_ylabel("$P(\omega)$",rotation=0)
    ax.yaxis.labelpad = 20
    ax.set_xlabel("$\omega$")
    ax.legend()

plt.tight_layout()

#fig.savefig("../plots/windowing_effects.eps")
#fig.savefig("../plots/beta_{:.2f}_spectrum.eps".format(B))
fig.savefig("../plots/multiple_spectra.eps".format(B))

plt.show()

