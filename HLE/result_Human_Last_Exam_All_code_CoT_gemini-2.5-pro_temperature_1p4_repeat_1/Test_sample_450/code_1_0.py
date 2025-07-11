import numpy as np

# Step 1: Define constants
epsilon_0 = 6.9e-3  # eV
T = 4200.0          # Kelvin
k_B = 8.61733e-5    # Boltzmann constant in eV/K
N_total_moles = 3.0 # Total number of particles in moles

# Thermal energy
kT = k_B * T

# Number of energy levels
num_levels = 5
n_values = np.arange(1, num_levels + 1)

# Step 2: Calculate energy levels and degeneracies
# Calculate energy increments epsilon_n
epsilons = epsilon_0 * (1 + n_values**n_values)

# Calculate total energy E_n (cumulative sum of epsilons)
energies = np.cumsum(epsilons)

# Calculate degeneracies g_n
degeneracies = 2 * n_values + 1

# Step 3: Compute the single-particle partition function (Z1)
boltzmann_factors = degeneracies * np.exp(-energies / kT)
Z1 = np.sum(boltzmann_factors)

print(f"The single-particle partition function Z1 is: {Z1:.4f}")

# Step 4: Determine the number of moles per level
# Calculate the probability of a particle being in each energy level
probabilities = boltzmann_factors / Z1

# Calculate the number of moles in each energy level
moles_per_level = N_total_moles * probabilities

# Step 5: Output the final result
print("The most probable number of particles in each energy level (in moles) is:")
# The problem asks to output each number in the final equation.
# The final equation is the tuple (m1, m2, m3, m4, m5)
final_tuple = tuple(moles_per_level)
print(final_tuple)

# Final answer in the required format
print(f"<<<{final_tuple}>>>")
