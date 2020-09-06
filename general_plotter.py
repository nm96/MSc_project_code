from simulation import Simulation
import sys
import pickle
import matplotlib.pyplot as plt

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)
    
ar = 3/2

fig = plt.figure(figsize = [6.4, ar*6.4])

ax1 = fig.add_axes([0.1, 0.1/ar, 0.8, 0.3]) # Dimensions as fractions of fig width/height.
ax2 = fig.add_axes([0.1, 0.3 + 0.2/ar, 0.3, 0.3])
ax2 = fig.add_axes([0.6, 0.3 + 0.2/ar, 0.3, 0.3])


plt.show()
