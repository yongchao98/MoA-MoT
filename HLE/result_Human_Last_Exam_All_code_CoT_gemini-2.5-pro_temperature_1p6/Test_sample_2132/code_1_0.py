# The product of the Favorskii rearrangement of 2-chlorocyclohexanone
# is ethyl cyclopentanecarboxylate.
# We will calculate its Böttcher Molecular Complexity.

# The formula for Böttcher Molecular Complexity (BMC) is:
# BMC = (N_atoms * N_bonds * N_cycles) / 2

# N_atoms: Number of non-hydrogen atoms.
# For ethyl cyclopentanecarboxylate (C8H14O2), this is 8 carbons + 2 oxygens.
N_atoms = 10

# N_bonds: Total number of bonds between non-hydrogen atoms (counting bond order).
# - 7 C-C single bonds
# - 2 C-O single bonds
# - 1 C=O double bond (counts as 2)
# Total N_bonds = 7 + 2 + 2 = 11
N_bonds = 11

# N_cycles: Number of rings in the molecule.
# There is one cyclopentane ring.
N_cycles = 1

# Calculate the Böttcher Molecular Complexity
bmc_value = (N_atoms * N_bonds * N_cycles) / 2

# Print the calculation with each number explicitly shown
print(f"The Böttcher Molecular Complexity is calculated as (N_atoms * N_bonds * N_cycles) / 2.")
print(f"Calculation: {N_atoms} * {N_bonds} * {N_cycles} / 2 = {bmc_value}")
