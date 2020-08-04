from simulation import Simulation
import pickle

s1 = Simulation()

# Set any non-default parameters here:

s1.B = 1.65

# Solve and save data:

s1.solve()

with open("latest.pkl", 'wb') as f:
    pickle.dump(s1,f,pickle.HIGHEST_PROTOCOL)

afn = "../data/B{:.2f}_Om{:.2f}_c{:.2f}_T{:.0e}.pkl"
afn = afn.format(s1.B, s1.Om, s1.c, s1.T)

with open(afn, 'wb') as f:
    pickle.dump(s1,f,pickle.HIGHEST_PROTOCOL)
