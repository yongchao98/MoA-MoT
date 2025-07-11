import math

# The problem asks for the limit of M(theta) as theta -> 0.
# Based on the analysis, this limit is a constant value determined by the geometry.
# The value is arccos(1/5).

# Step 1: The quantity we need to find is the limit of the supremum of the angle
# of incidence alpha.
# M = lim_{theta->0} sup(alpha)

# Step 2: The analysis shows that sup(alpha) is related to inf(gamma), where gamma is
# the angle between the trajectory and the side A, by sup(alpha) = pi/2 - inf(gamma).

# Step 3: Using the unfolding principle for billiards, the problem is reduced to
# finding the minimum angle between a line tangent to the unit circle and the line of side A.

# Step 4: In the limit as theta -> 0, side A becomes the line x=5.

# Step 5: The minimum angle, inf(gamma), occurs for a line from (5,0) tangent to the
# unit circle. In the right triangle formed by (0,0), (5,0), and the tangent point,
# we find that sin(inf(gamma)) = 1/5.

# Step 6: This gives inf(gamma) = arcsin(1/5).

# Step 7: Therefore, the final answer is M = pi/2 - arcsin(1/5) = arccos(1/5).

# Let's define the ratio for the calculation.
ratio = 1/5

# Calculate arccos(1/5) in radians.
result_rad = math.acos(ratio)

# The result is arccos(1/5)
print(f"The limit of M(theta) as theta goes to 0 is arccos(1/5).")
print(f"The equation for the final answer is M = arccos(1 / 5).")
print(f"Numerically, this value is approximately: {result_rad}")
