# The product of the Favorskii rearrangement of 2-chlorocyclohexanone is cyclopentanecarboxylic acid.
# The chemical formula is C6H10O2.

# The Böttcher Molecular Complexity is calculated as: (N_atoms * N_bonds) / N_parts

# 1. N_atoms: Total number of atoms in the molecule.
# For C6H10O2, N_atoms = 6 (C) + 10 (H) + 2 (O) = 18
n_atoms = 18

# 2. N_parts: Total number of non-hydrogen atoms.
# N_parts = 6 (C) + 2 (O) = 8
n_parts = 8

# 3. N_bonds: Total number of atomic connections (multiple bonds count as 1).
# For a single molecule with one ring, N_bonds = N_atoms.
# So, N_bonds = 18
n_bonds = 18

# Calculate the Böttcher Molecular Complexity Index (BMCI)
bmci = (n_atoms * n_bonds) / n_parts

# Print the final equation with the calculated values
print(f"The Böttcher Molecular Complexity is calculated as: (N_atoms * N_bonds) / N_parts")
print(f"For cyclopentanecarboxylic acid:")
print(f"({n_atoms} * {n_bonds}) / {n_parts} = {bmci}")
