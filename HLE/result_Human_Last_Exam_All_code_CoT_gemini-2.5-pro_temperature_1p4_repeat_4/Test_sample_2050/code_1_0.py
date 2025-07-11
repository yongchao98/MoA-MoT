import numpy as np

# System parameters
N = 7
num_photons = 4

# The ground state in the hard-core boson limit (U -> infinity) is found by
# filling the lowest single-particle energy levels of the equivalent
# tight-binding Hamiltonian for non-interacting fermions.

# The single-particle energies for a ring of N sites are given by:
# E_m = -2 * J * cos(2 * pi * m / N)
# where m are integer quantum numbers.

# To find the ground state for 4 particles, we fill the 4 states with the
# lowest energy. By inspecting the values of cos(2*pi*m/7), we find the
# lowest energy states correspond to modes m = 0, 1, -1, and 2.
occupied_modes = [0, 1, -1, 2]

# Calculate the energy contribution from each occupied mode (in units of J)
energies_J = [-2 * np.cos(2 * np.pi * m / N) for m in occupied_modes]

# Calculate the total hopping energy by summing the individual energies
total_hopping_energy_J = sum(energies_J)

# The total Hamiltonian is H = w * (total photons) + H_hopping
# So, E_ground = num_photons * w + E_hopping
# E_ground = 4*w + (E_0 + E_1 + E_-1 + E_2)

print("The Hamiltonian is H = H_const + H_hopping.")
print(f"For {num_photons} photons, the constant energy term is H_const = {num_photons}*ω.")
print("\nThe hopping energy is the sum of the energies of the 4 lowest single-particle states.")
print(f"These states correspond to the mode numbers m = {occupied_modes}.")
print("\nThe energy contribution from each mode (E_m = -2*J*cos(2πm/N)) is:")

# Output each number required for the final equation
print(f"E(m={occupied_modes[0]}) = {energies_J[0]:.4f} * J")
print(f"E(m={occupied_modes[1]}) = {energies_J[1]:.4f} * J")
print(f"E(m={occupied_modes[2]}) = {energies_J[2]:.4f} * J")
print(f"E(m={occupied_modes[3]}) = {energies_J[3]:.4f} * J")
print("-" * 30)
# Use a + sign for the sum to be explicit, even if the value is negative
print(f"Total hopping energy = {total_hopping_energy_J:.4f} * J")
print("\nThe final expression for the ground state energy E is:")

# Print the final equation with all components, as requested
print(f"E = {num_photons}*ω + ({energies_J[0]:.4f} + {energies_J[1]:.4f} + {energies_J[2]:.4f} + {energies_J[3]:.4f})*J")
print(f"E = {num_photons}*ω + ({total_hopping_energy_J:.4f})*J")
