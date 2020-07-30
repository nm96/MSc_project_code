import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib import rc
import sys
import pickle

# Set a nice font for the plots when on Linux
if sys.platform == 'linux':
    rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    rc('text', usetex=True)


inpf = sys.argv[1] # i.e. "spec_data.pkl"

with open(inpf, 'rb') as f:
    params, om, P = pickle.load(f)

# Note:
# model = (mdl.NHSommerfeld2,(Om,h,B))
# params = (eps,Om,m,c,k,model)

om_max = om[-1]
Om = params[1]
om_nat = (params[4]/params[2])**0.5
B = params[-1][-1][-1]
c = params[3]

fn = 0  # Initialize figure number for plotting
fn += 1; fig = plt.figure(fn,figsize=[12,6])

# Plot spectrum, label:
ax = fig.add_subplot(111)
ax.set_xlim([0,om_max])
ax.axvline(om_nat,ls='--',c='g',label=r"$\omega_{nat}$")
ax.axvline(Om,ls='--',c='r',label=r"$\Omega$")
ax.semilogy(om,P,c='k')
locmaj = matplotlib.ticker.LogLocator(base=100,numticks=30) 
ax.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=100,subs=(0.2,0.4,0.6,0.8),numticks=50)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.grid()
ax.set_title(r"""$\beta = {:.2f}$, $\Omega = {:.2f}$, $c = {:.2f}$""".format(B,Om,c))
ax.set_ylabel("$P(\omega)$",rotation=0)
ax.yaxis.labelpad = 20
ax.set_xlabel("$\omega$")
ax.legend()

plt.tight_layout()

#fig.savefig("../plots/sparse_peak_spectrum.eps".format(B))

plt.show()

