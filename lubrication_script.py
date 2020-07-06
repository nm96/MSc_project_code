import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfft
import matplotlib.pyplot as plt

from solutions import *
from models import *

fn = 0  # Initialize figure number for plotting

# Rename basic functions and constants for clarity
cos = np.cos
sin = np.sin
tanh = np.tanh
pi = np.pi

# Define parameter values:

eps = 0.0525 # Rotor eccentricity
Om = 4.1 # Driving frequency
m = 10 # Mass (per unit length)
c = 0.5 # Damping coefficient
k = 10 # Stiffness coefficient
h = 0.1 # Gap width
k_c = 50 # Stator stiffness parameter for models
alpha = 500 # for the Hua model
mu = 0.1 # Coefficient of friction

Om_nat = (k/m)**0.5 # Shaft natural frequency
