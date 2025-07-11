import math

# This script prints the final derived formula for the cutoff frequency.

# The Thevenin resistance R_th found from the ladder network analysis is:
# R_th = r * (1 + sqrt(3))

# The angular cutoff frequency omega_c is given by the formula:
# omega_c = 1 / (R_th * C)

print("The cutoff frequency (in radians per second) is given by the equation:")
print("omega_c = 1 / (R_th * C)")
print("\nSubstituting the derived value for the Thevenin resistance R_th = r * (1 + sqrt(3)), we get the final equation:")
print("omega_c = 1 / (r * (1 + sqrt(3)) * C)")
print("\nTo show each number in the final equation as requested:")
numerator = 1
num_in_factor = 1
num_in_sqrt = 3
print(f"Final Equation: omega_c = {numerator} / (r * C * ({num_in_factor} + sqrt({num_in_sqrt})))")

# Provide an approximate numerical factor for context
factor = 1 + math.sqrt(3)
print(f"\nThe numerical factor (1 + sqrt(3)) is approximately {factor:.4f}.")