# A regular icosahedron has 12 vertices and 30 edges.
# The great icosahedron shares these vertices, so its convex hull is a regular icosahedron.
# The surface of a regular icosahedron is already made of triangles.
# We can find the number of these triangular faces (F) using Euler's polyhedron formula: V - E + F = 2
# Where V is the number of vertices and E is the number of edges.

# Number of vertices of a regular icosahedron
V = 12
# Number of edges of a regular icosahedron
E = 30
# Euler's characteristic for a sphere
chi = 2

# We calculate the number of faces (F) using the formula: F = chi - V + E
F = chi - V + E

print(f"The outer-hull of a great icosahedron is a regular icosahedron.")
print(f"A regular icosahedron has {V} vertices and {E} edges.")
print("Using Euler's polyhedron formula (V - E + F = 2), we can find the number of faces (F).")
print("F = 2 - V + E")
print(f"F = {chi} - {V} + {E}")
print(f"The number of triangular faces is: {F}")