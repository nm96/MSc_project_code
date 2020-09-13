from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)

s1.phase_plot(d=2)

fig = plt.gcf()

#fig.suptitle("Phase-space plots for simulation at standard parameter values",y=0.9)

# Figure saving options:

#plt.savefig("../final_report/figures/standard_phase_plot.png")
#plt.savefig("../final_report/figures/standard_phase_plot.pdf")
plt.savefig("../final_report/figures/nonperiodic_phase_plot.pdf")

plt.show()
