# Step 1: Define the properties of the product molecule, cyclopentanecarboxylic acid.

# Number of non-hydrogen atoms (6 Carbon + 2 Oxygen)
num_heavy_atoms = 8

# Total number of bonds (connections) in the molecule
num_bonds = 18

# Number of rings
num_rings = 1

# Number of double bonds (the C=O in the carboxyl group)
num_double_bonds = 1

# Number of stereocenters (the ring carbon attached to the carboxyl group)
num_stereocenters = 1

# Step 2: Calculate the Böttcher Molecular Complexity using the defined formula.
# Complexity = N_heavy_atoms + N_bonds + N_rings + N_double_bonds + N_stereocenters
bottcher_complexity = num_heavy_atoms + num_bonds + num_rings + num_double_bonds + num_stereocenters

# Step 3: Print the final equation and the result.
print("The Böttcher Molecular Complexity is calculated as the sum of its components:")
print(f"Complexity = (Heavy Atoms) + (Bonds) + (Rings) + (Double Bonds) + (Stereocenters)")
print(f"Complexity = {num_heavy_atoms} + {num_bonds} + {num_rings} + {num_double_bonds} + {num_stereocenters}")
print(f"Total Böttcher Molecular Complexity: {bottcher_complexity}")