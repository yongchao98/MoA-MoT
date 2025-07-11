import math
from fractions import Fraction

# Step 1: Define initial parameters and calculate the area of the first circle (the first term, 'a').
w0 = Fraction(6)
h0 = Fraction(8)

# The radius of the first circle is one-third of the width.
r1 = w0 / 3

# The area of the first circle is the first term of the series. We represent it as a_val * pi.
a_val = r1**2

# Step 2: Calculate the scaling factor to find the common ratio.
# We find the intersection points of the first circle with the diagonals.
# The rectangle is centered at (0,0), so vertices are at (+-w0/2, +-h0/2).
# The equation of the diagonal is y = (h0/w0)*x. The circle is x^2 + y^2 = r1^2.
# Solving for the x-coordinate of the intersection:
# x_intersect = r1 * w0 / sqrt(w0^2 + h0^2) = 2 * 6 / sqrt(36+64) = 12 / 10 = 6/5.
x_intersect = Fraction(6, 5)

# The new smaller rectangles are created in the corners. The width of a new rectangle, w1, is
# found by taking the original half-width (w0/2) and subtracting the intersection's x-coordinate.
# w1 = w0/2 - x_intersect = 6/2 - 6/5 = 3 - 6/5 = 9/5.
# Or, more formally, w1 = (w0 - 2 * x_intersect) / 2
w1 = (w0 - 2 * x_intersect) / 2

# The linear scaling factor k is the ratio of the new width to the old width.
k = w1 / w0

# Step 3: Calculate the common ratio 'r' of the geometric series.
# At each step, we add 4 new circles, and their area scales by k^2.
# So, the common ratio of the areas added at each step is r = 4 * k^2.
r_val = 4 * k**2

# Step 4: Calculate the sum of the infinite geometric series and print the equation steps.
print("The total area of all circles is the sum of an infinite geometric series, S = a / (1 - r).")
print(f"The first term (area of the initial circle) is a = {a_val.numerator}*pi.")
print(f"The common ratio is r = {k.numerator**2*4}/{k.denominator**2} = {r_val.numerator}/{r_val.denominator}.")
print("\n--- Calculation ---")

# Substitute the values into the formula
print(f"S = ({a_val.numerator} * pi) / (1 - {r_val.numerator}/{r_val.denominator})")

# Calculate 1 - r
one_minus_r = 1 - r_val
print(f"S = ({a_val.numerator} * pi) / ({one_minus_r.numerator}/{one_minus_r.denominator})")

# Invert and multiply
final_num_unsimplified = a_val.numerator * one_minus_r.denominator
final_den_unsimplified = one_minus_r.numerator
print(f"S = ({final_num_unsimplified} * pi) / {final_den_unsimplified}")

# Simplify the final fraction by dividing by the greatest common divisor
common_divisor = math.gcd(final_num_unsimplified, final_den_unsimplified)
simplified_num = final_num_unsimplified // common_divisor
simplified_den = final_den_unsimplified // common_divisor
print(f"S = ({simplified_num}*pi) / {simplified_den}")