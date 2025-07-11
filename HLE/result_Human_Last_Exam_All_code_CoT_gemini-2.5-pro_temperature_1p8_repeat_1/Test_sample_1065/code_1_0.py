import math

# This script calculates the coordinates of the center of mass for the candy string.
# Based on the problem description, we model the string as an arc on a spherical
# quadrant of radius R. For symmetry, we assume the string follows the path that
# bisects the quadrant, from the top pole to the equator (φ = π/4).

# The formulas for the horizontal coordinates of the center of mass are derived as:
# X_cm = (√2 * R) / π
# Y_cm = (√2 * R) / π

# The problem asks for the "raw number" of these coordinates.
# We interpret this as the numerical coefficient of R in the expressions.
coefficient = math.sqrt(2) / math.pi

# The final output should be the horizontal coordinates (x and y) separated by a comma.
# In the equation for the center of mass, the numbers are the coefficient, the radius R, and pi.
# As requested, we provide the calculated numerical coefficient for each coordinate.
print(f"{coefficient},{coefficient}")