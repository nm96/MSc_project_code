import numpy as np
import solutions as sln
import pickle
import sys

inpf = sys.argv[1] # i.e. "sol_data.pkl"

with open(inpf, 'rb') as f:
    params, sol = pickle.load(f)

om_max = 10 # Maximum frequency to include in the spectral data

om, P = sln.PSD(sol)

P = P[om < om_max]
om = om[om < om_max]

spec_data = (params, om, P)

with open("spec_data.pkl", 'wb') as f:
    pickle.dump(spec_data,f,pickle.HIGHEST_PROTOCOL)
    
B = params[-1][-1][-1]
Om = params[1]
c = params[3]
T = sol.t[-1]

fname = "../data/B{:.2f}_Om{:.2f}_c{:.2f}_T{:.0e}_spec_data.pkl".format(B,Om,c,T)
with open(fname, 'wb') as f:
    pickle.dump(spec_data,f,pickle.HIGHEST_PROTOCOL)
