import math

# The derived expression for the cutoff angular frequency wc is:
# wc = 1 / ( (1 + sqrt(3)) * r * C )
# This script prints the formula with the numerical coefficients evaluated,
# as requested by the prompt.

# Define the constants from the derived formula
numerator = 1
constant_in_denominator_1 = 1
constant_in_denominator_2 = 3
symbolic_r = 'r'
symbolic_C = 'C'

# Calculate the numerical value of the coefficient for (r * C)
numerical_coefficient = 1 + math.sqrt(3)

print("The determined cutoff angular frequency (wc) at node a0 is:")
# Print the final equation with symbolic and numerical values
print(f"wc = {numerator} / (({constant_in_denominator_1} + sqrt({constant_in_denominator_2})) * {symbolic_r} * {symbolic_C})")
print("\nNumerically, this is approximately:")
print(f"wc = {numerator} / ({numerical_coefficient:.5f} * {symbolic_r} * {symbolic_C})")
