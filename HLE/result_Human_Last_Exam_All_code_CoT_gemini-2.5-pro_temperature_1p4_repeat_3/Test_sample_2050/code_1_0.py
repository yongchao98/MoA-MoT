import numpy as np

# System parameters
N = 7  # Number of cavities
M = 4  # Number of photons

# The ground state energy for the effective Hamiltonian is given by:
# E_ground = M * w + E_hopping
# where E_hopping is the sum of the M lowest single-particle energies.

# The single-particle energies are E_q = -2*J*cos(2*pi*q/N).
# The M=4 lowest energy levels for N=7 correspond to q = 0, 1, -1, 2.

# Calculate the values of the cosine terms
cos_q0 = np.cos(0)
cos_q1 = np.cos(2 * np.pi / N)
cos_q2 = np.cos(4 * np.pi / N)

# The sum of the cosine terms for the occupied levels
sum_of_cosines = cos_q0 + 2 * cos_q1 + cos_q2

# The coefficient of the J term in the energy
j_coefficient = -2 * sum_of_cosines

# Print the derivation of the ground state energy
print("The ground state energy is E_ground = 4*w + E_hopping.")
print("E_hopping is calculated by summing the energies of the M=4 lowest states.")
print(f"The occupied states are q = 0, 1, -1, 2.")
print("E_hopping = -2*J * (cos(2*pi*0/7) + cos(2*pi*1/7) + cos(2*pi*(-1)/7) + cos(2*pi*2/7))")
print("E_hopping = -2*J * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))")
print("\nCalculating the numerical values:")
print(f"cos(0) = {cos_q0:.5f}")
print(f"cos(2*pi/7) = {cos_q1:.5f}")
print(f"cos(4*pi/7) = {cos_q2:.5f}")
print("\nSubstituting the values into the equation for the J-dependent part:")
print(f"E_hopping = -2*J * ({cos_q0:.5f} + 2 * {cos_q1:.5f} + {cos_q2:.5f})")
print(f"E_hopping = -2*J * ({sum_of_cosines:.5f})")
print(f"E_hopping = {j_coefficient:.5f} * J")
print("\nCombining with the constant w term, the final ground state energy is:")
print(f"E_ground = 4*w + ({j_coefficient:.5f})*J")
print(f"E_ground = 4\u03C9 - {-j_coefficient:.5f}J")
