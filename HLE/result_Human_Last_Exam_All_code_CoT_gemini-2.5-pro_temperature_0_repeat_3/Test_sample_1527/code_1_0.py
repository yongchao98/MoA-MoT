import math

# Step 1 & 2: Define the geometric properties.
# The two lines are L1: y = x + 1 (or x - y + 1 = 0) and L2: y = -x + 5 (or x + y - 5 = 0).
# The lines are perpendicular and intersect at (2, 3).
# The angle bisectors are the lines x = 2 and y = 3.
# The center of any circle tangent to both lines must lie on one of these bisectors.

# Let circle C have center (cx, cy) and radius r.
# Let circle D have center (dx, dy) and radius 2.

# The radius R of a circle with center on a bisector is given by the distance formula to one of the lines.
# For a center (x, 3) on the bisector y=3, the radius R is |x - 2| / sqrt(2).
# For a center (2, y) on the bisector x=2, the radius R is |y - 3| / sqrt(2).

# Step 3 & 4: Analyze the configuration.
# As reasoned in the plan, the only valid configuration for two distinct tangent circles is
# for their centers to lie on different angle bisectors.
# Let's place center C on the bisector y = 3, and center D on the bisector x = 2.
# C = (cx, 3)
# D = (2, dy)

# Step 5: Set up the equations.
# For circle C: r = |cx - 2| / sqrt(2)  => (cx - 2)^2 = 2 * r^2
# For circle D: 2 = |dy - 3| / sqrt(2)  => (dy - 3)^2 = 2 * 2^2 = 8

# The distance between the centers C and D squared is:
# dist_sq = (cx - 2)^2 + (3 - dy)^2
# Since the circles are tangent, the distance between their centers is the sum of their radii, r + 2.
# dist_sq = (r + 2)^2

# Substituting the expressions from the radius formulas:
# (r + 2)^2 = (2 * r^2) + 8

# Step 6: Solve the equation for r^2.
# We need to solve the equation: r^2 + 4r + 4 = 2r^2 + 8
# This simplifies to: r^2 - 4r + 4 = 0
# This is a quadratic equation in the form of a perfect square.

# We can represent the coefficients of the final quadratic equation r^2 - 4r + 4 = 0
a = 1
b = -4
c = 4

# The equation is (r - 2)^2 = 0, which gives r = 2.
r = 2

# The problem asks for the value of r^2.
r_squared = r**2

print("The final equation to solve for the radius r is:")
print(f"{a}r^2 + ({b})r + {c} = 0")
print("This can be factored as (r - 2)^2 = 0.")
print(f"Solving for r gives r = {r}.")
print(f"The value of r^2 is {r}^2, which is:")
print(r_squared)
