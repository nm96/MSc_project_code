# A script for integrating the Jeffcott equations using scipy's solve_ivp.

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Rename basic functions and constants for clarity:
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi
