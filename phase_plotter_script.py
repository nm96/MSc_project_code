from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)

s1.phase_plot(d=2)

# Save commands
#plt.savefig("../report/B1_65phaseplot.png",dpi=400)
#plt.savefig("../report/B1_5phaseplot.png",dpi=400)
#plt.savefig("../report/Om2_7phaseplot.png",dpi=400)
#plt.savefig("../report/Om0_3phaseplot.png",dpi=400)

plt.show()
