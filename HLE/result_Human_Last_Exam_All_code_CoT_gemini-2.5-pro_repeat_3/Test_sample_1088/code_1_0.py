# Define the atom for calculation
atom = "Carbon"

# Define the components of the 6-311G** basis set for a second-row heavy atom
# Core shell (1s) is represented by a contraction of 6 primitives
core_primitives = 6

# Valence shells (2s, 2p) are triple-split (311)
# Each valence orbital is represented by 3+1+1 primitives
valence_primitives_per_orbital = 3 + 1 + 1

# Carbon has one 2s orbital and three 2p orbitals in its valence shell
num_s_valence_orbitals = 1
num_p_valence_orbitals = 3

# The first '*' in ** adds d-polarization functions to heavy atoms (5 functions)
d_polarization_primitives = 5

# The second '*' adds p-polarization functions to Hydrogen, so it's 0 for Carbon.
p_polarization_primitives = 0

# Calculate the total for each component
s_valence_total = num_s_valence_orbitals * valence_primitives_per_orbital
p_valence_total = num_p_valence_orbitals * valence_primitives_per_orbital

# Calculate the final total number of primitive Gaussians
total_primitives = core_primitives + s_valence_total + p_valence_total + d_polarization_primitives

# Print the breakdown of the calculation
print(f"Calculation for a {atom} atom with the 6-311G** basis set:")
print(f"Core primitives (1s): {core_primitives}")
print(f"Valence s-primitives (2s): {s_valence_total}")
print(f"Valence p-primitives (2p): {p_valence_total}")
print(f"Polarization d-primitives (*): {d_polarization_primitives}")
print("\nFinal Equation:")
print(f"{core_primitives} + {s_valence_total} + {p_valence_total} + {d_polarization_primitives} = {total_primitives}")
