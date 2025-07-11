import math

# The problem asks for the ratio Area(D)/Area(S) in its exact form.
# The side length 's' cancels out in the final ratio. Our derivation below assumes s=1 for simplicity.
# The final fraction is independent of s.

# The reasoning for the final expression is as follows:
# 1. Total Surface Area, Area(S) = 6 * s^2.
# 2. The region D consists of points on the surface with distance <= sqrt(2)*s from a vertex P.
# 3. On the 3 faces adjacent to P, all points are within this distance. Area = 3 * s^2.
# 4. On each of the 3 faces not adjacent to P, the area of D is s^2 * (pi/3 + (sqrt(3)-3)/2).
# 5. Total Area(D) = 3*s^2 + 3 * s^2 * (pi/3 + (sqrt(3)-3)/2)
#                   = s^2 * (3 + pi + 3*sqrt(3)/2 - 9/2)
#                   = s^2 * (2*pi + 3*sqrt(3) - 3) / 2.
# 6. The ratio Area(D)/Area(S) is [s^2*(2*pi+3*sqrt(3)-3)/2] / [6*s^2].
#    This simplifies to (2*pi + 3*sqrt(3) - 3) / 12.

# The code below presents the coefficients of this final exact expression.
# The expression is of the form: (A*pi + B*sqrt(3) + C) / D

numerator_pi_coefficient = 2
numerator_sqrt3_coefficient = 3
numerator_constant_term = -3
denominator = 12

print("The final answer is a fraction based on an exact geometric calculation.")
print("The final equation is in the form: (A*pi + B*sqrt(3) + C) / D")
print("\nThe values of the integer coefficients and the denominator are:")
print(f"A (coefficient of pi): {numerator_pi_coefficient}")
print(f"B (coefficient of sqrt(3)): {numerator_sqrt3_coefficient}")
print(f"C (constant term): {numerator_constant_term}")
print(f"D (denominator): {denominator}")

print("\nThus, the final ratio is:")
# We display the final formula with the calculated numbers
print(f"({numerator_pi_coefficient}*pi + {numerator_sqrt3_coefficient}*sqrt(3) {numerator_constant_term}) / {denominator}")
