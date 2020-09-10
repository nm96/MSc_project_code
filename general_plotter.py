from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf, opf = sys.argv[1], sys.argv[2]

with open(inpf,'rb') as f:
    s1 = pickle.load(f)
    
ar = 5/4

fig = plt.figure(figsize = [6.4, ar*6.4])

ax1 = fig.add_axes([0.2, 0.1/ar, 0.7, 0.20]) # Dimensions as fractions of fig width/height.
ax2 = fig.add_axes([0.2, 0.22 + 0.2/ar, 0.7, 0.20])
ax3 = fig.add_axes([0.2, 0.5 + 0.2/ar, 0.3, 0.3/ar])
ax4 = fig.add_axes([0.6, 0.5 + 0.2/ar, 0.3, 0.3/ar])

om_max = 5  # Remember to change this to an appropriate value!

s1.phase_plot(axes=(ax3,ax4))
s1.psd_plot(ax = ax1,om_max = om_max)
s1.Rpsd_plot(ax = ax2, om_max = om_max)

ax1.set_title(r"PSD plot")
ax2.set_title(r"RPSD plot")

fig.text(0.41, 0.92, r"Phase space trajectories",fontsize=12)

ts = r"""Simulation results for $\Omega = {:.3f}$, $\beta={:.3f}$"""
ts = ts.format(s1.Om, s1.B)
fig.suptitle(ts,fontsize=14)

plt.savefig(opf)
#plt.savefig("../final_report/figures/test_plot.png")
#plt.savefig("../final_report/figures/test_plot.pdf")

plt.show()
