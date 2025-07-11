# Step 1: Define the properties of the product molecule, methyl cyclopentanecarboxylate (C7H12O2).

# Number of atoms of each element
n_carbon = 7
n_hydrogen = 12
n_oxygen = 2

# Total number of atoms (N_atoms)
n_atoms = n_carbon + n_hydrogen + n_oxygen

# Number of covalent bonds (N_bonds), calculated from atomic valences (C=4, H=1, O=2)
# N_bonds = (sum of valences) / 2
n_bonds = (n_carbon * 4 + n_hydrogen * 1 + n_oxygen * 2) // 2

# Number of cycles (N_cycles) in the molecule
n_cycles = 1

# Number of disconnected parts (N_parts) of the molecule
n_parts = 1

# Step 2: Calculate the Böttcher Molecular Complexity using the formula.
# BCI = (N_atoms * N_bonds * N_cycles) / N_parts
boettcher_complexity = (n_atoms * n_bonds * n_cycles) / n_parts

# Step 3: Print the result and the equation used.
print("The product of the Favorskii rearrangement of 2-chlorocyclohexanone is methyl cyclopentanecarboxylate.")
print("Its Böttcher Molecular Complexity is calculated as follows:")
print(f"Number of atoms (N_atoms) = {n_atoms}")
print(f"Number of bonds (N_bonds) = {n_bonds}")
print(f"Number of cycles (N_cycles) = {n_cycles}")
print(f"Number of parts (N_parts) = {n_parts}")
print("\nFinal Calculation:")
# The final line prints the equation with the calculated values
print(f"{n_atoms} * {n_bonds} * {n_cycles} / {n_parts} = {int(boettcher_complexity)}")