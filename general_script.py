from simulation import Simulation
import matplotlib.pyplot as plt

s1 = Simulation()

s1.Om = 2.7; s1.c = 0.2; s1.B = 2.2
s1.T = 2**10

s1.solve()
s1.transform()
s1.find_peaks()

for om in s1.peaks:
    print("\t $$ \t & \t {:.3f} \t \\\\".format(om))

s1.psd_plot()

plt.show()
