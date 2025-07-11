# The problem is to find the maximum number of grid squares the perimeter of a triangle
# with sides 18, 18, and 18*sqrt(2) can cross without touching any lattice points.

# This can be modeled by finding three integer vectors, v1, v2, v3, that represent
# the sides of the triangle on a grid. These vectors must satisfy:
# 1. Their Euclidean lengths are 18, 18, and 18*sqrt(2).
# 2. They form a closed loop: v1 + v2 + v3 = (0, 0).
# The number of squares crossed, k, is the sum of their Manhattan lengths (|x| + |y|).

# By analysis, the integer vectors that maximize the Manhattan length for the given
# Euclidean lengths are those aligned with the axes or diagonals.
# For a side of length 18, the vector is (18, 0) or (0, 18). Manhattan length = 18.
# For a side of length 18*sqrt(2), the vector is (18, 18). Manhattan length = 36.

# We choose a set of vectors that satisfies the conditions:
v1 = (18, 0)
v2 = (-18, 18)
v3 = (0, -18)

# Number of squares crossed by the first leg (length 18)
k1 = abs(v1[0]) + abs(v1[1])

# Number of squares crossed by the hypotenuse (length 18*sqrt(2))
k2 = abs(v2[0]) + abs(v2[1])

# Number of squares crossed by the second leg (length 18)
k3 = abs(v3[0]) + abs(v3[1])

# Total number of squares crossed
k_total = k1 + k2 + k3

print(f"The number of squares crossed by the first leg is: {k1}")
print(f"The number of squares crossed by the second leg is: {k3}")
print(f"The number of squares crossed by the hypotenuse is: {k2}")
print(f"The total number of squares crossed is the sum of these values.")
print(f"{k1} + {k3} + {k2} = {k_total}")
print(f"The largest number k is: {k_total}")
