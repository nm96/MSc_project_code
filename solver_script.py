from simulation import Simulation
import pickle

s1 = Simulation()

# Set any non-default parameters here:

s1.Om = 1
s1.B = 20.45
s1.c = 50
s1.k = 0.00
s1.T = 2**10

# The following effectively removes m from the equations without affecting
# anything:
#s1.m *= 0.1
#s1.c *= 0.1
#s1.k *= 0.1
#s1.B *= 0.1

# Solve and save data:

s1.solve()

with open("latest.pkl", 'wb') as f:
    pickle.dump(s1,f,pickle.HIGHEST_PROTOCOL)

afn = "../data/B{:.3f}_Om{:.2f}_c{:.2f}_T{:.0e}.pkl"
afn = afn.format(s1.B, s1.Om, s1.c, s1.T)

if s1.T > 2**10:
    with open(afn, 'wb') as f:
        pickle.dump(s1,f,pickle.HIGHEST_PROTOCOL)
