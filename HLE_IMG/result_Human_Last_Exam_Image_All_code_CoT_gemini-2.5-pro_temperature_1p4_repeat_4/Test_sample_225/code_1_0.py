# Define atomic weights (g/mol)
atomic_weight_C = 12.011
atomic_weight_H = 1.008
atomic_weight_O = 15.999

# Define the number of atoms of each element in Compound A (C19H16O7)
num_C = 19
num_H = 16
num_O = 7

# Calculate the molecular weight
molecular_weight = (num_C * atomic_weight_C) + (num_H * atomic_weight_H) + (num_O * atomic_weight_O)

# Print the final equation and the result
print(f"The molecular formula of compound A is C{num_C}H{num_H}O{num_O}.")
print("The calculation for the molecular weight is:")
print(f"({num_C} * {atomic_weight_C}) + ({num_H} * {atomic_weight_H}) + ({num_O} * {atomic_weight_O}) = {molecular_weight:.3f}")
print(f"The molecular weight of compound A is {molecular_weight:.3f} g/mol.")