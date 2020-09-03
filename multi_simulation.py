import simulation
import numpy as np
import pickle

c_set = np.linspace(0.2,0.6,21)
periods = []

for c in c_set:
    s = simulation.Simulation()
    s.B = 1.65
    s.c = c
    s.T = 2**12
    s.solve()
    s.find_period()
    periods.append(s.Tp)

print(c_set,periods)
