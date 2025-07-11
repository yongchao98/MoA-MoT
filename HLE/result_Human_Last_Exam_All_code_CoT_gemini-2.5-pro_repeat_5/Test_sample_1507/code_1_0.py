# The outer-hull of a great icosahedron is a regular icosahedron.
# We need to find the number of its faces.

# A regular icosahedron has 12 vertices.
V = 12

# 5 edges meet at each vertex of a regular icosahedron.
k = 5

# The total number of edges (E) can be calculated using the number of vertices.
# Since each edge is shared by two vertices, we have: E = (V * k) / 2
E = (V * k) // 2

# Using Euler's formula for polyhedra (V - E + F = 2), we can find the number of faces (F).
# F = 2 - V + E
F = 2 - V + E

print("The number of faces (F) of the regular icosahedron (the outer-hull) is calculated using Euler's formula: F = 2 - V + E")
print(f"Number of Vertices (V) = {V}")
print(f"Number of Edges (E) = ({V} * 5) / 2 = {E}")
print("The final equation is:")
print(f"F = 2 - {V} + {E}")
print(f"F = {F}")
print(f"\nThe minimal number of triangles needed is {F}.")