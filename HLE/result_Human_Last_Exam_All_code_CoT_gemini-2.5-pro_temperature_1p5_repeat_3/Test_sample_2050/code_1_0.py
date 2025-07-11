import numpy as np

# Description of the final energy equation and its components
# The total ground state energy 'E' is given by the formula:
# E = c_omega * omega + c_J * J
# where c_omega is a simple integer and c_J is a more complex coefficient.

# System parameters
N = 7  # Number of cavities
M = 4  # Number of photons

# The coefficient of omega is simply the number of photons.
c_omega = M

# The coefficient of J is derived from the ground state energy of the hopping term.
# Based on the analytical derivation, this coefficient can be expressed as:
# c_J = A * (B + C * cos(D * pi / F))
# We will define these numerical components.
A = -2
B = 1
C = 2
D = 2
F = N

# Print the final equation with all its numerical components explicitly shown.
print("The ground state energy E is given by the equation:")
print(f"E = {c_omega} * omega + ({A}) * ( {B} + {C} * cos({D}*pi/{F}) ) * J")

# For practical purposes, we can also compute the numerical value of the coefficient of J.
c_J_numeric = A * (B + C * np.cos(D * np.pi / F))
print("\nThe numerical value for the coefficient of J is approximately:")
print(f"{c_J_numeric:.5f}")
