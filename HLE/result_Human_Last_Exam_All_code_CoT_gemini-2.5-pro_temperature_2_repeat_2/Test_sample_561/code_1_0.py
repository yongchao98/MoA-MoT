import math

# The Minkowski-Bouligand dimension (D) for a self-affine fractal
# can be found by analyzing its scaling properties.
#
# Based on the problem description, we found:
# - The number of self-similar pieces is N = 5.
# - The scaling factor in the x-direction is s_x = 1/14.
# - The scaling factor in the y-direction is s_y = 9/14.
#
# The growth rate of the number of covering boxes scales with N * (s_y / s_x),
# which is 5 * ( (9/14)/(1/14) ) = 5 * 9 = 45.
# The shrinking rate of the box size scales with 1/s_x = 14.
#
# This leads to the formula D = log(45) / log(14).

numerator_in_log = 45
denominator_in_log = 14

# Calculate the dimension.
# The base of the logarithm (e.g., 10, e) does not matter
# as log_b(x) / log_b(y) is constant for any base b.
dimension = math.log(numerator_in_log) / math.log(denominator_in_log)

# Print the final formula and the calculated result.
print("The formula for the dimension (D) is:")
print(f"D = log({numerator_in_log}) / log({denominator_in_log})")
print("\nThe calculated dimension of the black keys is:")
print(dimension)