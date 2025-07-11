import math

# The problem of finding the Minkowski–Bouligand dimension for the self-affine
# set of black keys can be solved using the box-counting method.

# Based on the derivation:
# N = 5 (the number of black keys to iterate on)
# s_x = 1/14 (the scaling factor for the width)
# s_y = 9/14 (the scaling factor for the height)

# The dimension D can be calculated as log(N_effective) / log(1/s_effective).
# Through a detailed box-counting analysis, the dimension D is given by the formula:
# D = log(45) / log(14)

numerator_val = 45
denominator_val = 14

# Calculate the dimension
dimension = math.log(numerator_val) / math.log(denominator_val)

# Output the equation and the result
print("The Minkowski–Bouligand dimension (D) is calculated using the formula:")
print(f"D = log({numerator_val}) / log({denominator_val})")
print(f"The calculated dimension is: {dimension}")