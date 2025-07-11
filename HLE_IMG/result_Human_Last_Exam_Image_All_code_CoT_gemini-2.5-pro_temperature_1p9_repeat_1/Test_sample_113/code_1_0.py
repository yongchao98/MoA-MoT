import math
from fractions import Fraction

# The problem asks for the limit of the total area of all circles as n -> infinity.
# This can be modeled as the sum of an infinite geometric series.
# Total Area S = A_gen1 + A_gen2 + A_gen3 + ...
# where A_gen_k is the total area of all circles introduced at step k.

# Step 1: Calculate the area of the first circle (the first term of the series, a).
w0 = 6
h0 = 8

# The radius of the first circle is one-third of the width.
r1_val = Fraction(w0, 3)

# The area of a circle is pi * r^2. We will calculate the coefficient of pi.
# This coefficient is the first term, 'a', of our series.
a = r1_val**2

print("The calculation for the total area S involves summing a geometric series: S = a / (1 - r).")
print("First, we find the first term 'a', which is the area of the initial circle.")
print(f"Initial rectangle width = {w0}")
print(f"Radius of the first circle, r1 = 1/3 * {w0} = {r1_val}")
print(f"Area of the first circle = pi * ({r1_val})^2 = {a}*pi")
print(f"So, the first term of the series is a = {a}*pi.\n")

# Step 2: Determine the common ratio 'r' of the series.
# This involves finding the scaling factor from one generation of circles to the next.

# At each step, 4 new smaller figures (and circles) are created.
num_new_circles = 4

# To find the scaling factor, we determine the dimensions of the new rectangles.
# The original rectangle's half-width and half-height are:
half_w0 = Fraction(w0, 2)
half_h0 = Fraction(h0, 2)

# The circle (radius r1) intersects the diagonals. We find the intersection point coordinates.
# Let the center be (0,0). Top-right vertex is at (3,4). Diagonal line is y = (4/3)x.
# Circle equation is x^2 + y^2 = r1^2 = 4.
# Substituting y gives: x^2 + (16/9)x^2 = 4 => (25/9)x^2 = 4 => x = 6/5.
intersect_x = Fraction(6, 5)

# The new rectangles at the corners are formed in the space between the intersection points
# and the original vertices. The width of a new corner rectangle is:
w1 = half_w0 - intersect_x

# The linear scaling factor for dimensions is the ratio of the new width to the old width.
linear_scale_factor = w1 / w0

# The area of a new circle scales by the square of the linear factor.
area_scale_factor = linear_scale_factor**2

# The common ratio 'r' is the total area of the new generation divided by the total area of the previous one.
# This is the number of new circles multiplied by their individual area scaling factor.
r = num_new_circles * area_scale_factor

print("Next, we find the common ratio 'r'.")
print(f"The linear scaling factor from one generation of rectangles to the next is {linear_scale_factor}.")
print(f"The area of a single circle thus scales by ({linear_scale_factor})^2 = {area_scale_factor}.")
print(f"Since {num_new_circles} new circles are created in each step, the total area added at each step is scaled by a factor of r = {num_new_circles} * {area_scale_factor} = {r}.")
print(f"So, the common ratio of the series is r = {r}.\n")


# Step 3: Sum the infinite geometric series using the formula S = a / (1 - r).
denominator = 1 - r
total_area_coeff = a / denominator

print("Finally, we calculate the sum of the infinite series.")
print("The final equation is S = a / (1 - r)")
print(f"Plugging in the numbers for the coefficients of pi: S_coeff = {a} / (1 - {r})")
print(f"S_coeff = {a} / ({denominator})")
print(f"S_coeff = {total_area_coeff}")
print(f"\nTherefore, the limit of the total area is {total_area_coeff.numerator}*pi / {total_area_coeff.denominator}.")
<<<25*pi/4>>>