import math

# Number of self-similar pieces (black keys)
N = 5

# The horizontal scaling ratio is 1/14.
# The vertical scaling ratio is 9/14.
# The dimension D is given by the formula log(N * r_y/r_x) / log(1/r_x)
# which simplifies to log(N * 9) / log(14)

numerator_arg = 45  # This comes from N * (r_y/r_x) = 5 * ( (9/14) / (1/14) ) = 5 * 9
denominator_arg = 14 # This comes from 1/r_x = 1 / (1/14)

dimension = math.log(numerator_arg) / math.log(denominator_arg)

print(f"The fractal dimension D is calculated by the formula: D = log({N} * 9) / log({denominator_arg})")
print(f"This simplifies to the equation: D = log({numerator_arg}) / log({denominator_arg})")
print(f"The calculated Minkowski-Bouligand dimension is: {dimension}")