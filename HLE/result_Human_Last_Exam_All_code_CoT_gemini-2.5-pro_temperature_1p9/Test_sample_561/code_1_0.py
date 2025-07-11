import math

# The Minkowski-Bouligand dimension for a self-affine fractal of this type
# can be calculated from the number of copies (N) and the scaling ratios (r_x, r_y).
# The formula is D = log(N * r_y / r_x) / log(1 / r_x)

# Based on the problem description:
# N = 5 (five black keys)
# r_x = 1/14 (scaling in the x-direction)
# r_y = 9/14 (scaling in the y-direction)

# We can simplify the terms in the formula.
# Numerator term: N * r_y / r_x = 5 * (9/14) / (1/14) = 5 * 9 = 45
numerator_base = 45

# Denominator term: 1 / r_x = 1 / (1/14) = 14
denominator_base = 14

# Calculate the dimension D = log(45) / log(14)
dimension = math.log(numerator_base) / math.log(denominator_base)

# Output the final equation and the result.
print("The formula for the dimension (D) is derived from the construction parameters.")
print(f"The final equation for the dimension is D = log({numerator_base}) / log({denominator_base})")
print("\nThis evaluates to:")
print(f"D = {math.log(numerator_base)} / {math.log(denominator_base)}")
print(f"D = {dimension}")
