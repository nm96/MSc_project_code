import simulation
import numpy as np
import pickle
import time
import datetime

t0 = time.time()

#periods = []

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

#B_set = np.linspace(1.2,1.3,41)
#for B in B_set:
#    s = simulation.Simulation()
#    s.c = 0.4
#    s.B = B
#    s.T = 2**12
#    s.solve()
#    s.find_period()
#    periods.append(s.Tp)
#
#print("Note c = 0.4 rather than the default 0.5")
#for k in range(len(B_set)):
#    print("B = {:.3f} --> Tp = {:.3f}, om1 = {:.3f}".format(B_set[k],
#        periods[k], 2*np.pi/periods[k]))

def period_at(B,c):
    s = simulation.Simulation()
    s.B = B
    s.c = c
    s.T = 2**12
    s.solve()
    s.find_period()
    return s.Tp

N = 11
M = 6

B_range = (1.2,1.3)
c_range = (0.4,0.7)

B_set = np.linspace(*B_range,N)
c_set = np.linspace(*c_range,M)

periods = np.zeros((N,M))

for n in range(N):
    for m in range(M):
        periods[n,m] = period_at(B_set[n],c_set[m])

print(periods)

run = {'B_set':B_set, "c_set":c_set, "periods":periods}

fn = "../data/multi_parameter_runs/Bcperiods.pkl"
with open(fn,'wb') as f:
    pickle.dump(run, f, pickle.HIGHEST_PROTOCOL)

tf = time.time()
print("Total time taken = ",datetime.timedelta(seconds = tf-t0))
