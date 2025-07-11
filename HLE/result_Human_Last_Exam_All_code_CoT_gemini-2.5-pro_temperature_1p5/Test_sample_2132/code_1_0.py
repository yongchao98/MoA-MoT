# The product of the Favorskii rearrangement of 2-chlorocyclohexanone
# is cyclopentanecarboxylic acid.

# Chemical Formula: C6H10O2

# Step 1: Calculate the total number of atoms (N_atoms).
# Carbons (C): 6
# Hydrogens (H): 10
# Oxygens (O): 2
num_atoms = 6 + 10 + 2

# Step 2: Calculate the total number of bonds (N_bonds) in the molecular graph.
# A bond, whether single, double, or triple, counts as one connection.
# C-C bonds: 6 (5 in the ring, 1 connecting the -COOH group)
# C-H bonds: 10
# C=O bond: 1
# C-O bond: 1
# O-H bond: 1
num_bonds = 6 + 10 + 1 + 1 + 1

# Step 3: Calculate the Böttcher Molecular Complexity.
# For a single molecule, Complexity = N_atoms * N_bonds.
bottcher_complexity = num_atoms * num_bonds

# Step 4: Print the final result and the equation.
print("The Böttcher Molecular Complexity for cyclopentanecarboxylic acid is calculated as:")
print(f"Number of Atoms * Number of Bonds = Complexity")
print(f"{num_atoms} * {num_bonds} = {bottcher_complexity}")