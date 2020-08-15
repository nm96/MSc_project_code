from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)

s1.psd_plot()

# Save commands
#plt.savefig("../report/B1_65phaseplot.png",dpi=400)
#plt.savefig("../report/Om2_7psdplot.png",dpi=400)
plt.savefig("../report/Om2_7apsdplot.png",dpi=400)
#plt.savefig("../report/Om0_3psdplot.png",dpi=400)

plt.show()
