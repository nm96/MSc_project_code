from simulation import Simulation
import matplotlib.pyplot as plt

s1 = Simulation()

s1.Om = 2.7; s1.c = 0.2; s1.B = 2.2
s1.T = 2**14

s1.first_return_plot()

plt.show()
