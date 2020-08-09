from simulation import Simulation
import matplotlib.pyplot as plt

s1 = Simulation()

s1.B = 1.65
s1.T = 2**12

s1.first_return_plot()

plt.show()
