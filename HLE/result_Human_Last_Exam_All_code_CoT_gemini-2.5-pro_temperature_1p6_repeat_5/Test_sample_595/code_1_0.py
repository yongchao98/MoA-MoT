import math

# We are tasked with finding the maximum number of grid squares 'k' that the perimeter of an
# isosceles right triangle (sides 18, 18, 18*sqrt(2)) can pass through without touching any
# lattice points.

# The strategy is to position the triangle optimally on the coordinate plane. We can model this by
# defining the vectors for the two legs of length 18. Let these vectors be (u, v) and (-v, u),
# which are perpendicular and have a magnitude of sqrt(u^2+v^2). We need this magnitude to be 18,
# so u^2 + v^2 = 18^2 = 324.

# To maximize the number of grid squares crossed, we need to maximize the number of grid lines
# crossed. By placing a vertex at a position like (epsilon, epsilon) where epsilon is a very
# small positive number, we can derive a formula for the number of crossed squares, k,
# in terms of U = floor(u) and V = floor(v).

# The formula, assuming U and V are positive, is k = 4 * max(U, V) + 2 * min(U, V) + 2.

# We must find non-negative integers U and V that can be realized, meaning the circle u^2+v^2=324
# must intersect the square [U, U+1) x [V, V+1). The conditions for this are:
# 1) U^2 + V^2 <= 324
# 2) (U+1)^2 + (V+1)^2 >= 324

# By searching through possible integer pairs (U,V) from (0,0) to (18,18), we find that the
# pair that maximizes the formula for k while satisfying the conditions is (16, 8) or (8, 16).
# We choose U=16 and V=8, where U is the larger value.

U = 16
V = 8

# Calculate the terms for the final equation
term1 = 4 * U
term2 = 2 * V
term3 = 2
max_k = term1 + term2 + term3

print(f"The optimal integer pair found is U={U}, V={V}.")
print("Using the formula k = 4*max(U,V) + 2*min(U,V) + 2, we can calculate the maximum number of squares.")
print("The final equation is:")
print(f"{4} * {U} + {2} * {V} + {2} = {term1} + {term2} + {term3} = {max_k}")