import numpy as np
import solutions as slns
import pickle

om_max = 10 # Maximum frequency to include in the spectral data

with open("sol_data.pkl", 'rb') as f:
    params, sol = pickle.load(f)

om, P = slns.PSD(sol)

P = P[om < om_max]
om = om[om < om_max]

spec_data = (params, om, P)

with open("spec_data.pkl", 'wb') as f:
    pickle.dump(spec_data,f,pickle.HIGHEST_PROTOCOL)
