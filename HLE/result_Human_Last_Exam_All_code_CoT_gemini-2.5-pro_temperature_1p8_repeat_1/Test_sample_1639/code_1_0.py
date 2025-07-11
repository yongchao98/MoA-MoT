# Masses of the most stable isotopes (monoisotopic masses)
mass_C = 12.000000
mass_H = 1.007825
mass_O = 15.994915

# The final modification on cysteine is determined to be a fragment
# with the chemical formula C5H10O4.
num_C = 5
num_H = 10
num_O = 4

# Calculate the total mass of the modification
modification_mass = (num_C * mass_C) + (num_H * mass_H) + (num_O * mass_O)

# Print the calculation steps
print(f"The mass of the modification is calculated from its formula: C{num_C}H{num_H}O{num_O}")
print(f"Mass = ({num_C} * Carbon_mass) + ({num_H} * Hydrogen_mass) + ({num_O} * Oxygen_mass)")
print(f"Mass = ({num_C} * {mass_C}) + ({num_H} * {mass_H}) + ({num_O} * {mass_O})")
print(f"Mass = {num_C * mass_C} + {num_H * mass_H:.4f} + {num_O * mass_O:.4f}")
print(f"Final Mass (x) = {modification_mass:.4f}")

# The question asks for the number x, which is the integer mass.
# Based on the calculation, the mass is ~134 Da.
print(f"\nThe number x is the integer value of the calculated mass.")
print(f"x = {round(modification_mass)}")
