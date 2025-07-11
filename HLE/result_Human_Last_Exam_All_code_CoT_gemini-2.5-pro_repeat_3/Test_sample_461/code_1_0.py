# Pterodactyladiene has the chemical formula C30H16.
# This script calculates its molecular weight.

# Define the atomic weights (in atomic mass units, u)
atomic_weight_C = 12.011
atomic_weight_H = 1.008

# Define the number of atoms for each element in the molecule
num_C = 30
num_H = 16

# Calculate the total mass from each element
total_weight_C = num_C * atomic_weight_C
total_weight_H = num_H * atomic_weight_H

# Calculate the total molecular weight
molecular_weight = total_weight_C + total_weight_H

# Print the name of the molecule
print("Molecule: Pterodactyladiene (C30H16)")
print("Calculation of Molecular Weight:")

# Print the equation with all the numbers used
print(f"({num_C} * {atomic_weight_C}) + ({num_H} * {atomic_weight_H}) = {molecular_weight:.3f} u")