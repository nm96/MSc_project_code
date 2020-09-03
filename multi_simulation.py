import simulation
import numpy as np
import pickle
import time
import datetime

t0 = time.time()

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

B_set = np.linspace(1.2,1.4,21)
for B in B_set:
    s = simulation.Simulation()
    s.B = B
    s.T = 2**12
    s.solve()
    s.find_period()
    periods.append(s.Tp)

for k in range(len(B_set)):
    print("B = {:.2f} --> Tp = {:.3f}, om1 = {:.3f}".format(B_set[k],
        periods[k], 2*np.pi/periods[k]))

tf = time.time()
print("Total time taken = ",datetime.timedelta(seconds = tf-t0))
