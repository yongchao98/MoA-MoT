#
# This script calculates the molecular weight of the final product, Compound 4.
#

# Chemical formula of 2,2'-benzoyldibenzoic acid is C15H10O5
num_C = 15
num_H = 10
num_O = 5

# Standard atomic weights
atomic_weight_C = 12.011  # g/mol
atomic_weight_H = 1.008   # g/mol
atomic_weight_O = 15.999  # g/mol

# Calculate the contribution of each element
mass_C = num_C * atomic_weight_C
mass_H = num_H * atomic_weight_H
mass_O = num_O * atomic_weight_O

# Calculate the total molecular weight
total_molecular_weight = mass_C + mass_H + mass_O

# Print the final equation and the result
print("The final product is 2,2'-benzoyldibenzoic acid (C15H10O5).")
print("Its molecular weight is calculated as follows:")
print(f"({num_C} * {atomic_weight_C}) + ({num_H} * {atomic_weight_H}) + ({num_O} * {atomic_weight_O}) = {total_molecular_weight:.3f} g/mol")
