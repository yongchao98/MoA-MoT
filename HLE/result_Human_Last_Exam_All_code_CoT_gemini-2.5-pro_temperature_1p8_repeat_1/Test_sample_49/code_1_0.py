# The f-vector of a polytope is a sequence of integers (f_0, f_1, f_2, ..., f_{d-1})
# where f_i is the number of i-dimensional faces.

# We are given a non-simplicial 4-polytope with 6 vertices.
# This means d=4 and f_0=6.
f_0 = 6

# The polytope has the maximal number of 2-faces (f_2).
# The Upper Bound Theorem for polytopes states that the maximum number of faces
# is achieved by the neighborly polytopes.
# For 4 dimensions and 6 vertices, there are two such polytopes:
# 1. The simplicial cyclic polytope C(6,4).
# 2. A unique non-simplicial neighborly polytope.
# Both have the same number of vertices, edges, and 2-faces.
# The f-vector of the boundary complex for both is (f_0, f_1, f_2) = (6, 15, 20).
f_1 = 15
f_2 = 20

# We can find the number of 3-faces (facets), f_3, using the Euler-Poincaré formula for 4-polytopes:
# f_0 - f_1 + f_2 - f_3 = 0
# Rearranging the formula to solve for f_3: f_3 = f_0 - f_1 + f_2
f_3 = f_0 - f_1 + f_2

# Print the components of the f-vector
print(f"The number of vertices is f_0 = {f_0}")
print(f"The number of edges is f_1 = {f_1}")
print(f"The number of 2-faces is f_2 = {f_2}")
print(f"The number of 3-faces (facets) is f_3 = {f_3}")
print("-" * 20)
print("This is derived from the Euler-Poincaré formula: f0 - f1 + f2 - f3 = 0")
print(f"The final equation is: {f_0} - {f_1} + {f_2} - {f_3} = {f_0 - f_1 + f_2 - f_3}")
print("-" * 20)
print(f"The f-vector of the polytope is ({f_0}, {f_1}, {f_2}, {f_3}).")
