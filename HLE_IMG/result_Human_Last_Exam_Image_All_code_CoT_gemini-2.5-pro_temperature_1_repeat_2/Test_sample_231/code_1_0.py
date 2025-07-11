# The final compound C is 1,3,6,8-tetrahydroxy-10,10-dimethyl-10H-pyrano[3,2-b]xanthen-5-ylium.
# The molecular formula of the cation is [C21H14O6]+.
# This script calculates its molecular mass using integer atomic masses.

# Atomic masses
mass_C = 12
mass_H = 1
mass_O = 16

# Number of atoms in the cation of compound C
num_C = 21
num_H = 14
num_O = 6

# Calculate the total mass
total_mass = (num_C * mass_C) + (num_H * mass_H) + (num_O * mass_O)

# Print the final equation with each number, as requested.
print("The molecular formula of the cation of compound C is C21H14O6.")
print("The equation to calculate the molecular mass is:")
print(f"({num_C} * {mass_C}) + ({num_H} * {mass_H}) + ({num_O} * {mass_O}) = {total_mass}")
print(f"The molecular mass of the cation of compound C is {total_mass} g/mol.")