from fractions import Fraction

# This script calculates the probability based on the analytical solution explained above.

# The expected area of a triangle formed by 3 random points in a unit square.
# This is a well-known result in geometric probability.
expected_area = Fraction(11, 144)

# The probability that the convex hull of 4 random points is a triangle (P_T).
# This is 4 times the expected area of the triangle, as any of the 4 points can be the one inside.
P_T = 4 * expected_area

# The probability 'p' that one point falls in the circumcircle of the other three
# is related to P_T by the formula: p = (2 - P_T) / 4.
p_numerator = 2 - P_T
p_denominator = 4
final_probability = p_numerator / p_denominator

# Print the final equation with all the numbers used in the calculation.
print("The probability 'p' is calculated by the formula: (2 - P_T) / 4")
print(f"First, we find P_T = 4 * (11/144) = {P_T.numerator}/{P_T.denominator}")
print("\nThen, we calculate the final probability:")
# The final line below shows each number in the final equation as requested.
print(f"p = (2 - {P_T.numerator}/{P_T.denominator}) / {p_denominator} = {final_probability.numerator}/{final_probability.denominator}")

<<<61/144>>>