import math

# The Minkowski–Bouligand dimension D for this self-affine set can be calculated
# using a box-counting argument.
# D = log(N_effective) / log(1/r_effective)
# Based on our analysis, the effective number of copies is N_effective = 45,
# and the effective scaling factor is r_effective = 1/14.

numerator = 45
denominator = 14

# Calculate the dimension
dimension = math.log(numerator) / math.log(denominator)

# Output the equation and the result
print(f"The dimension D is given by the formula: D = log({numerator}) / log({denominator})")
print(f"D = {math.log(numerator)} / {math.log(denominator)}")
print(f"D ≈ {dimension}")
