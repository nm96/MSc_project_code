import numpy as np
import solutions as sln
import pickle
import sys
import matplotlib.pyplot as plt

# Get input data file from command line
inpf = sys.argv[1] # i.e. "sol_data.pkl"

with open(inpf, 'rb') as f:
    params, sol = pickle.load(f)

Om = params[1]
h = params[-1][-1][1]

fig, ax = plt.subplots()

ax.plot(*sln.rotsolxy(sol,Om))
ax.plot(h*np.cos(np.linspace(0,2*np.pi,1000)),h*np.sin(np.linspace(0,2*np.pi,1000)),c='r') 
ax.set_aspect('equal')

#fig.savefig("../plots/phase_space_plot.eps")

plt.show()
