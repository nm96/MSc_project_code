from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)

s1.psd_plot(R=True)

ax = plt.gca()


s1.find_peaks()

n = 1
for om in s1.peaks:
    print("\t ${}\\tilde{{\omega}}_1$ \t & \t ${:.3f}$ \t \\\\".format(n,om))
    n += 1


# Save commands
#plt.savefig("../report/B1_65Rpsdplot.png",dpi=400)
#plt.savefig("../report/Om2_7psdplot.png",dpi=400)
#plt.savefig("../report/Om2_7apsdplot.png",dpi=400)
#plt.savefig("../report/Om0_3psdplot.png",dpi=400)

plt.show()

