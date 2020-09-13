from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)

s1.Fpsd_plot(om_max=15)


# Save commands

plt.savefig("../final_report/figures/standard_Fpsd_plot.png",dpi=200)
#plt.savefig("../final_report/figures/standard_Fpsd_plot.pdf")

#plt.savefig("../report/B1_65psdplot.png",dpi=400)
#plt.savefig("../report/B1_5psdplot.png",dpi=400)
#plt.savefig("../report/Om2_7psdplot.png",dpi=400)
#plt.savefig("../report/Om2_7apsdplot.png",dpi=400)
#plt.savefig("../report/Om0_3psdplot.png",dpi=400)

plt.show()

