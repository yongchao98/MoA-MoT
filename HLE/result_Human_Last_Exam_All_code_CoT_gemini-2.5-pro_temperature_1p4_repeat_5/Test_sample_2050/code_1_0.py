import numpy as np

# This script calculates the ground state energy for the given Hamiltonian
# in the limit of infinitely strong on-site repulsion (U -> infinity).

# Parameters of the system
N = 7  # Number of cavities
M = 4  # Number of photons

# In the U -> infinity limit, the system is equivalent to non-interacting fermions.
# The total energy is E_g = E_omega + E_J.
# E_omega = M * omega = 4*omega.
# E_J is the sum of the M=4 lowest single-particle energies of the hopping Hamiltonian.
# The single-particle energies are epsilon_m = -2*J*cos(2*pi*m/N).
# For M=4 on N=7 sites, the lowest energy states correspond to m = 0, +1, -1, and +2.
# So, E_J = -2*J * [cos(2*pi*0/7) + cos(2*pi*1/7) + cos(2*pi*(-1)/7) + cos(2*pi*2/7)]
# Which simplifies to:
# E_J = -2*J * [cos(0) + 2*cos(2*pi/7) + cos(4*pi/7)]

# Numerically evaluate the terms
c0 = np.cos(2 * np.pi * 0 / N)  # cos(0)
c1 = np.cos(2 * np.pi * 1 / N)  # cos(2*pi/7)
c2 = np.cos(2 * np.pi * 2 / N)  # cos(4*pi/7)

# Calculate the value inside the bracket
bracket_val = c0 + 2*c1 + c2

# The total coefficient for the J term
j_coeff = -2 * bracket_val

# Print the step-by-step calculation and the final result
print("The ground state energy E_g is the sum of the on-site energy and the hopping energy:")
print(f"E_g = {M}*omega + E_J")
print("The hopping energy E_J is a multiple of J. Its coefficient is calculated as follows:")
print(f"Coefficient = -2 * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))")
print(f"            = -2 * ({c0:.4f} + 2*({c1:.4f}) + ({c2:.4f}))")
print(f"            = -2 * ({bracket_val:.4f})")
print(f"            = {j_coeff:.4f}")
print("\nTherefore, the final equation for the ground state energy is:")
print(f"E_g = {M}*omega + {j_coeff:.4f}*J")
