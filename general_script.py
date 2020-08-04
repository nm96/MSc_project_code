from simulation import Simulation

s1 = Simulation()

s1.B = 1.65
s1.T = 2**10

s1.solve()
s1.phase_plot(1)

s1.transform()
s1.psd_plot(2)

plt.show()
