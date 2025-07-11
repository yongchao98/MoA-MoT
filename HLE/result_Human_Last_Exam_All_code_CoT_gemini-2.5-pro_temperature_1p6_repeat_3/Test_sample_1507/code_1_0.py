# The great icosahedron shares its 12 vertices with the regular icosahedron.
# Therefore, its convex hull (the "outer-hull") is a regular icosahedron.
# We want to find the number of triangular faces on the surface of this regular icosahedron.

# We can use Euler's formula for polyhedra: V - E + F = 2
# V = number of vertices
# E = number of edges
# F = number of faces

# For a regular icosahedron:
V = 12  # Number of vertices
E = 30  # Number of edges

# We need to solve for F (the number of faces/triangles).
# F = 2 - V + E
F = 2 - V + E

# Print the explanation and the result of the calculation.
print("The outer-hull of the great icosahedron is a regular icosahedron.")
print("Using Euler's formula (V - E + F = 2) to find the number of faces (F):")
print(f"Given V = {V} and E = {E}")
print(f"The equation is: {V} - {E} + F = 2")
print(f"Solving for F: F = 2 - {V} + {E}")
print(f"The number of triangles is: {F}")
