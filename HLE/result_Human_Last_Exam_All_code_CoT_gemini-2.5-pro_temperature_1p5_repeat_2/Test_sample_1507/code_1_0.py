# The great icosahedron is a non-convex polyhedron. Its convex hull,
# which represents the outer visible shape, is a regular icosahedron.
# We want to find the number of triangular faces on this regular icosahedron.
# We can use Euler's formula for polyhedra: V - E + F = 2
# where V is the number of vertices, E is the number of edges, and F is the number of faces.

# For a regular icosahedron:
# Number of vertices (V)
V = 12
# Number of edges (E)
E = 30

# We rearrange Euler's formula to solve for F: F = 2 - V + E
F = 2 - V + E

# Print the explanation and the calculation
print("The outer-hull of the great icosahedron is a regular icosahedron.")
print("To find the number of faces (triangles), we can use Euler's formula: V - E + F = 2.")
print("For a regular icosahedron, V = 12 and E = 30.")
print("Solving for F gives the equation: F = 2 - V + E")
print(f"Plugging in the values: F = 2 - {V} + {E}")
print(f"The calculation results in: F = {F}")
print(f"\nThus, the minimal number of triangles needed is {F}.")
