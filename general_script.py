from simulation import Simulation
import matplotlib.pyplot as plt
import sys
import pickle
import numpy as np

inpf = sys.argv[1] # Get input file from commmand line

with open(inpf,'rb') as f:
    s1 = pickle.load(f)

s1.find_period()

print(s1.Tp)
print(2*np.pi/s1.Tp)
