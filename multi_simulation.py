import simulation
import numpy as np
import pickle

periods = []

#c_set = np.linspace(0.2,0.6,21)
#for c in c_set:
#    s = simulation.Simulation()
#    s.B = 1.65
#    s.c = c
#    s.T = 2**12
#    s.solve()
#    s.find_period()
#    periods.append(s.Tp)
#
#for k in range(len(c_set)):
#    print("c = {:.2f} --> Tp = {:.4f}".format(c_set[k], periods[k]))

B_set = np.linspace(1.3,1.7,21)
for B in B_set:
    s = simulation.Simulation()
    s.B = B
    s.T = 2**12
    s.solve()
    s.find_period
    periods.append(s.Tp)

for k in range(len(c_set)):
    print("B = {:.2f} --> Tp = {:.4f}".format(B_set[k], periods[k]))
