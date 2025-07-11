import math

# Define the parameters of the system
N = 7  # Number of cavities in the ring
num_photons = 4  # Total number of photons

# The problem simplifies to finding the ground state energy of 4 non-interacting
# fermions on a 7-site ring. The total energy is E = 4*omega + E_hopping.
# We need to calculate E_hopping.

# The single-particle energies are E_m = -2J * cos(2*pi*m/N).
# The 4 particles fill the lowest 4 available energy levels, which correspond
# to momenta m = 0, 1, -1, and 2.

# The total hopping energy is the sum of these four single-particle energies:
# E_hopping = E(m=0) + E(m=1) + E(m=-1) + E(m=2)
# E_hopping = -2J*cos(0) - 2J*cos(2*pi/7) - 2J*cos(-2*pi/7) - 2J*cos(4*pi/7)
# Using cos(-x) = cos(x), this can be written as:
# E_hopping = -2J * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7))

# Calculate the numerical values of the cosine terms
cos_0 = math.cos(0)
cos_2pi_over_7 = math.cos(2 * math.pi / N)
cos_4pi_over_7 = math.cos(4 * math.pi / N)

# Calculate the final coefficient for the J term
J_coefficient = -2 * (cos_0 + 2 * cos_2pi_over_7 + cos_4pi_over_7)

# Print the final equation with all numbers shown
print("The ground state energy E is given by the equation:")
print(f"E = {num_photons}*omega + [ -2 * (cos(0) + 2*cos(2*pi/7) + cos(4*pi/7)) ] * J")

print("\nSubstituting the numerical values for the cosine terms, we have:")
# The format string below shows each number used in the calculation of the coefficient for J.
print(f"E = {num_photons} * omega + [ -2 * ( {cos_0:.6f} + 2*({cos_2pi_over_7:.6f}) + ({cos_4pi_over_7:.6f}) ) ] * J")

print("\nSimplifying the expression gives the final equation:")
print(f"E = {num_photons} * omega + {J_coefficient:.6f} * J")