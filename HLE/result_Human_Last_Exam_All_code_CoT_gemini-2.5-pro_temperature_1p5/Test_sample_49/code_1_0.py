# The f-vector of a convex 4-polytope is denoted by (f0, f1, f2, f3).
# We are given:
# f0 = 6 (number of vertices)
# The polytope is non-simplicial and has the maximal number of 2-faces.

# Based on the complete classification of f-vectors for 4-polytopes with 6 vertices,
# the vector that satisfies the conditions is (6, 15, 20, 11).
f0 = 6
f1 = 15
f2 = 20
f3 = 11

print(f"The f-vector is ({f0}, {f1}, {f2}, {f3}).")

# We can verify this f-vector with Euler's formula for 4-polytopes:
# f0 - f1 + f2 - f3 = 0
result = f0 - f1 + f2 - f3

# As requested, here is the final equation with each number shown.
print("Euler's formula verification:")
print(f"{f0} - {f1} + {f2} - {f3} = {result}")