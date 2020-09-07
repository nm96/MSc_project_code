from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)
    
ar = 4/3

fig = plt.figure(figsize = [6.4, ar*6.4])

ax1 = fig.add_axes([0.2, 0.1/ar, 0.7, 0.2]) # Dimensions as fractions of fig width/height.
ax2 = fig.add_axes([0.2, 0.23 + 0.2/ar, 0.7, 0.2])
ax3 = fig.add_axes([0.2, 0.55 + 0.2/ar, 0.3, 0.3/ar])
ax4 = fig.add_axes([0.6, 0.55 + 0.2/ar, 0.3, 0.3/ar])

s1.phase_plot(axes=(ax3,ax4))
s1.psd_plot(ax = ax1)
s1.psd_plot(R=True, ax = ax2)


ts = r"""Simulation Results With Parameters
$\beta={:.3f}$, $\Omega={:.2f}$, $c={:.2f}$"""
ts = ts.format(s1.B,s1.Om,s1.c)
fig.suptitle(ts)

#plt.savefig("../final_report/figures/test_plot.png")

plt.show()
