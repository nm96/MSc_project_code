"""Module containing functions for solving the Jeffcott rotordynamic equations,
processing and plotting the resulting solutions. Specific nonlinear models to
be incorporated in the Jeffcott equations are in the the 'models' module.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import fft
import matplotlib.pyplot as plt

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi


