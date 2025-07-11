import math

# Based on the analysis, we assume the center of the second circle (D) is at the
# intersection of the two tangent lines, and the tangency between circles is external.

# 1. Find the intersection of y = x + 1 and y = -x + 5.
# x + 1 = -x + 5  => 2x = 4 => x = 2
# y = 2 + 1 = 3
# Center of the second circle D is (2, 3).
# Radius of the second circle is given as 2.
R_d = 2

# 2. The center of the first circle, C, lies on an angle bisector (x=2 or y=3).
# Let's assume C = (2, cy). The distance from C to D is |cy - 3|.
# The radius r of the first circle is its distance to a tangent line, e.g., x + y - 5 = 0.
# r = |2 + cy - 5| / sqrt(2) = |cy - 3| / sqrt(2).
# This gives a key relationship: |cy - 3| = r * sqrt(2).

# 3. Use the external tangency condition for the two circles.
# The distance between centers C and D is the sum of their radii.
# distance(C, D) = r + R_d
# |cy - 3| = r + 2

# 4. Substitute and solve for r.
# From step 2 and 3, we have:
# r * sqrt(2) = r + 2
# r * sqrt(2) - r = 2
# r * (sqrt(2) - 1) = 2
# r = 2 / (sqrt(2) - 1)

# To simplify, we rationalize the denominator:
# r = 2 * (sqrt(2) + 1) / ((sqrt(2) - 1) * (sqrt(2) + 1))
# r = 2 * (sqrt(2) + 1) / (2 - 1)
# r = 2 * sqrt(2) + 2

# 5. Calculate r^2.
# r^2 = (2*sqrt(2) + 2)^2
# r^2 = (2*sqrt(2))^2 + 2*(2*sqrt(2))*2 + 2^2
# r^2 = 8 + 8*sqrt(2) + 4
# r^2 = 12 + 8*sqrt(2)

# The final equation for r^2 is r^2 = 12 + 8 * sqrt(2).
# The numbers in this equation are 12, 8, and 2.
a = 12
b = 8
c = 2

# Print the final equation as requested.
print(f"The final equation is: r^2 = {a} + {b}*sqrt({c})")

# Calculate and print the numerical value of r^2.
r_squared_value = a + b * math.sqrt(c)
print(f"The numerical value of r^2 is: {r_squared_value}")