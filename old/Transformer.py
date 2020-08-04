import numpy as np
import solutions as sln
import pickle
import sys

# Get input data file from command line
inpf = sys.argv[1] # i.e. "sol_data.pkl"

with open(inpf, 'rb') as f:
    params, sol = pickle.load(f)

om_max = 10 # Maximum frequency to include in the spectral data

# Use pre-defined function to calculate transform (could move function to here
# for simplicity?)
om, P = sln.PSD(sol)

# Cut spectrum at max frequency
P = P[om < om_max]
om = om[om < om_max]

# Bundle into a tuple 
spec_data = (params, om, P)

# Save spectrum to local binary file
with open("spec_data.pkl", 'wb') as f:
    pickle.dump(spec_data,f,pickle.HIGHEST_PROTOCOL)
    
# Get parameters in order to label the archive file
# (this step in particular would be a lot less messy with specialised objects?)
B = params[-1][-1][-1]
Om = params[1]
c = params[3]
T = sol.t[-1]

# Save external labelled file
fname = "../data/B{:.2f}_Om{:.2f}_c{:.2f}_T{:.0e}_specl_data.pkl".format(B,Om,c,T)
with open(fname, 'wb') as f:
    pickle.dump(spec_data,f,pickle.HIGHEST_PROTOCOL)
