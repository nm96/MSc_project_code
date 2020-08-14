from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)

#s1.psd_plot(R=True,om_max=20)
s1.psd_plot()

plt.savefig("../report/B0_5psdplot.png",dpi=400)

plt.show()
