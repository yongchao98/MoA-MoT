import numpy as np

def linf_norm(v):
    """
    Calculates the L-infinity norm of a vector.
    ||v||_inf = max(|v_i|)
    """
    return np.max(np.abs(v))

# Step 1: Establish an upper bound for the length.
# Let V be a normed real vector space. Let a and b be two points on the unit sphere.
# By definition of the unit sphere, we have ||a|| = 1 and ||b|| = 1.
# The length of the line segment connecting a and b is given by ||a - b||.
# Using the triangle inequality, which is a property of any norm:
# ||a - b|| = ||a + (-b)|| <= ||a|| + ||-b||
# Since ||-b|| = |-1| * ||b|| = ||b||, we have:
# ||a - b|| <= ||a|| + ||b|| = 1 + 1 = 2.
# This shows that the length of the segment can be at most 2.

# Step 2: Show that this upper bound of 2 is achievable.
# We will use the space R^2 with the L-infinity norm.
# The unit ball in this space, {x in R^2 : ||x||_inf <= 1}, is the square
# [-1, 1] x [-1, 1], which is a metrically convex set.
# The unit sphere is the boundary of this square.

# Let's choose two points, a and b, on the unit sphere.
# We choose the endpoints of one of the sides of the square.
a = np.array([1, 1])
b = np.array([1, -1])

# Verify that a and b are on the unit sphere by checking their norms.
norm_a = linf_norm(a)
norm_b = linf_norm(b)

print(f"Let V be the vector space R^2 with the L-infinity norm.")
print(f"Consider the points a = {a} and b = {b}.")
print(f"The norm of a is ||{a}||_inf = {norm_a}.")
print(f"The norm of b is ||{b}||_inf = {norm_b}.")
print("Since their norms are 1, points a and b are on the unit sphere.")

# We also need to confirm that the entire line segment connecting a and b lies on the unit sphere.
# A point on the segment is x(t) = (1-t)a + t*b for t in [0, 1].
# x(t) = (1-t)*[1, 1] + t*[1, -1] = [1, 1-2t].
# The norm is ||x(t)||_inf = max(|1|, |1-2t|).
# Since t is in [0, 1], 1-2t is in [-1, 1], so |1-2t| <= 1.
# Thus, ||x(t)||_inf = 1 for all t in [0, 1], confirming the segment is on the unit sphere.

print("\nThe line segment connecting a and b lies entirely on the unit sphere.")

# Finally, calculate the length of this segment.
diff_vector = a - b
length = linf_norm(diff_vector)

print("\nThe length of the segment is the norm of the difference between a and b.")
print(f"Length = ||a - b|| = ||{a} - {b}|| = ||{diff_vector}||_inf = {length}")

print("\nSince the maximum possible length is 2, and we have found a case where the length is 2, the largest possible length is 2.")