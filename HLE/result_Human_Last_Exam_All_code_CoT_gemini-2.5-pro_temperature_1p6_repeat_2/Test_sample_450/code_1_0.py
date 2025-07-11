import numpy as np

# Step 1: Define Physical Constants
epsilon_0_eV = 6.9e-3  # eV
T_K = 4200.0  # Kelvin
k_B_eV_K = 8.617333e-5  # Boltzmann constant in eV/K
total_moles = 3.0 # Total moles of particles

# Energy level indices
n_values = np.arange(1, 6)

# Step 2: Calculate Energy Level Properties
# Degeneracy of each level
g_n = 2 * n_values + 1
# Energy offsets for each level
epsilon_n = epsilon_0_eV * (1 + n_values**n_values)

# Step 3: Calculate Total Energy Levels
# E_n is the cumulative sum of epsilon_n
E_n = np.cumsum(epsilon_n)

# Thermal energy
kBT = k_B_eV_K * T_K

# Step 4: Calculate the Single-Particle Partition Function (Z_1)
# Boltzmann factor for each level, multiplied by degeneracy
terms = g_n * np.exp(-E_n / kBT)
Z1 = np.sum(terms)

# Step 5: Calculate the Number of Moles per Level
# The probability of a particle being in a given state
probabilities = terms / Z1
# The most probable number of moles in each energy level
moles_in_each_level = total_moles * probabilities

# Step 6: Output the Result
print(f"The single-particle partition function Z_1 is: {Z1:.4f}")
print("\nThe most probable number of moles in each energy level (E_1, E_2, E_3, E_4, E_5) is:")
# The problem asks for the equation, which we interpret as showing the final values.
print(f"({moles_in_each_level[0]:.4f}, {moles_in_each_level[1]:.4f}, {moles_in_each_level[2]:.4f}, {moles_in_each_level[3]:.4f}, {moles_in_each_level[4]:.4f})")

final_answer = tuple(moles_in_each_level)
# print(f"<<<{final_answer}>>>")