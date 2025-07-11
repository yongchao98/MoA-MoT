import math

# Step 1: Define the radii of the two circles.
# The second circle (center D) has a known radius.
r_d = 2
# The first circle (center C) has an unknown radius r.

# Step 2: Explain the geometric setup and key relationships.
# The lines y = x + 1 and y = -x + 5 intersect at I(2, 3).
# The lines are perpendicular, and their angle bisectors are x=2 and y=3.
# The center of any circle tangent to both lines lies on one of these bisectors.
# The distance from a circle's center to the intersection point I is given by:
# dist(Center, I) = radius * sqrt(2).
# So, dist(C, I) = r * sqrt(2) and dist(D, I) = r_d * sqrt(2).

# We assume the most symmetric configuration where the centers C and D lie on
# different bisectors (e.g., C on x=2 and D on y=3).
# This forms a right-angled triangle with vertices C, D, and I.

# Step 3: Use the Pythagorean theorem to find the squared distance between centers.
# dist(C, D)^2 = dist(C, I)^2 + dist(D, I)^2
# dist(C, D)^2 = (r * sqrt(2))^2 + (r_d * sqrt(2))^2
# dist(C, D)^2 = 2 * r^2 + 2 * r_d^2

# Step 4: Use the tangency condition between the two circles.
# For the circles to be tangent externally, the distance between their centers
# must be the sum of their radii.
# dist(C, D) = r + r_d
# Squaring this gives: dist(C, D)^2 = (r + r_d)^2

# Step 5: Equate the two expressions for dist(C, D)^2 and solve for r.
# 2 * r^2 + 2 * r_d^2 = (r + r_d)^2
# Substitute r_d = 2:
# 2 * r^2 + 2 * (2^2) = (r + 2)^2
# 2 * r^2 + 8 = r^2 + 4r + 4
# Rearranging gives the final quadratic equation for r:
# r^2 - 4r + 4 = 0

# Step 6: Print the coefficients of the final equation and solve for r^2.
a = 1
b = -4
c = 4
print(f"The final equation for r is: {a}*r^2 + ({b})*r + {c} = 0")
print(f"This is a perfect square: (r - 2)^2 = 0")
# The equation (r - 2)^2 = 0 has a single solution for r.
r = 2
print(f"The radius of the first circle is r = {r}")

# Calculate the final answer, r^2.
r_squared = r**2
print(f"The value of r^2 is: {r_squared}")

<<<4>>>