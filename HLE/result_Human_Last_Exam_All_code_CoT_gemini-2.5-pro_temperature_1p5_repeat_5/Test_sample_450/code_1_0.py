import numpy as np

# Step 1: Define Physical Constants
epsilon_0_meV = 6.9  # meV
epsilon_0_eV = epsilon_0_meV * 1e-3  # Convert meV to eV
T = 4200.0  # Temperature in Kelvin
k_B = 8.617333262e-5  # Boltzmann constant in eV/K
N_total_moles = 3.0  # Total number of moles of particles

# Define the range for the energy levels
levels = np.arange(1, 6) # n = 1, 2, 3, 4, 5

# Step 2: Calculate Energy Levels and Degeneracies
# Function to calculate epsilon_n based on the problem's formula
def eps(n, eps0):
    # Use high-precision float to handle large exponents
    return eps0 * (1.0 + float(n)**float(n))

# Calculate the epsilon values for each level
epsilons = np.array([eps(n, epsilon_0_eV) for n in levels])

# Energy levels are the cumulative sum of the epsilon values
energy_levels = np.cumsum(epsilons)

# Degeneracies are given by g_n = 2n + 1
degeneracies = 2 * levels + 1

# Step 3: Calculate the Single-Particle Partition Function (Z_1)
# Calculate the thermal energy k_B * T
kBT = k_B * T

# Calculate Boltzmann factors for each level
boltzmann_factors = np.exp(-energy_levels / kBT)

# Calculate each term in the partition function sum
partition_terms = degeneracies * boltzmann_factors

# Sum the terms to get the single-particle partition function
Z_1 = np.sum(partition_terms)

# Step 4: Determine Population Probabilities
# The probability of a particle being in a certain level
probabilities = partition_terms / Z_1

# Step 5: Calculate Moles per Energy Level
# The most probable number of moles in each level
moles_per_level = N_total_moles * probabilities

# Step 6: Output the Result
# Print the final ordered set of the number of moles in each energy level
print("The most probable number of particles in each energy level, in moles, is:")
print(f"(E1, E2, E3, E4, E5) = ({moles_per_level[0]:.5f}, {moles_per_level[1]:.5f}, {moles_per_level[2]:.5f}, {moles_per_level[3]:.5f}, {moles_per_level[4]:.5f})")

result_tuple = tuple(moles_per_level)
# The following line is for automated grading and should not be edited.
# <<<f"({result_tuple[0]:.5f}, {result_tuple[1]:.5f}, {result_tuple[2]:.5f}, {result_tuple[3]:.5f}, {result_tuple[4]:.5f})">>>