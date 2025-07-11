# The outer-hull (convex hull) of a great icosahedron is a regular icosahedron.
# We need to find the number of triangular faces on a regular icosahedron.
# We can use Euler's formula for polyhedra: V - E + F = 2
# where V = vertices, E = edges, and F = faces.

# A regular icosahedron has 12 vertices.
V = 12

# A regular icosahedron has 30 edges.
E = 30

# We rearrange Euler's formula to solve for F:
# F = 2 - V + E
F = 2 - V + E

print("The outer-hull of the great icosahedron is a regular icosahedron.")
print(f"A regular icosahedron has {V} vertices (V) and {E} edges (E).")
print("Using Euler's formula (F = 2 - V + E) to find the number of faces:")
print(f"Number of Triangles = 2 - {V} + {E}")
print(f"The minimal number of triangles needed is: {F}")