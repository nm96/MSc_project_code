from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt
import matplotlib.ticker

inpf1 = sys.argv[1] # Get input file from commmand line
inpf2 = sys.argv[2] # Get input file from commmand line

with open(inpf1,'rb') as f:
    s1 = pickle.load(f)
    
with open(inpf2,'rb') as f:
    s2 = pickle.load(f)

s1.transform(R=True)
s2.transform(R=True)

fig, ax = plt.subplots(figsize=[6,4])

ax.set_xlim([0,10])
ax.semilogy(s1.om,s1.P,label=r"$\beta = {:.2f}$".format(s1.B))
ax.semilogy(s2.om,s2.P,c='red',label=r"$\beta = {:.2f}$".format(s2.B))
locmaj = matplotlib.ticker.LogLocator(base=100,numticks=30) 
ax.yaxis.set_major_locator(locmaj)
locmin = matplotlib.ticker.LogLocator(base=100,subs=(0.2,0.4,0.6,0.8),
        numticks=50)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.grid()
ax.set_ylabel("$P(\omega)$",rotation=0)
ax.yaxis.labelpad = 20
ax.set_xlabel("$\omega$")
ax.legend()
plt.tight_layout()

# Save commands
#plt.savefig("../report/B1_65psdplot.png",dpi=400)
#plt.savefig("../report/B1_5psdplot.png",dpi=400)
#plt.savefig("../report/Om2_7psdplot.png",dpi=400)
#plt.savefig("../report/Om2_7apsdplot.png",dpi=400)
#plt.savefig("../report/Om0_3psdplot.png",dpi=400)
#plt.savefig("../report/Cpsdplot.png",dpi=400)

plt.show()

