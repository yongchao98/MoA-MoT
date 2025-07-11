# The outer-hull of the great icosahedron is a regular icosahedron.
# We need to find the number of its faces (F), which are all triangles.
# We use Euler's formula for polyhedra: V - E + F = 2
# For a regular icosahedron:
# Number of vertices (V)
V = 12
# Number of edges (E)
E = 30

# Solve for F: F = 2 - V + E
F = 2 - V + E

print("The outer-hull is a regular icosahedron. We find its number of faces (F) using Euler's formula: V - E + F = 2.")
print(f"For an icosahedron, V = {V} and E = {E}.")
print("The equation with the final number of faces is:")
print(f"{V} - {E} + {F} = 2")
print(f"\nThus, the minimal number of triangles needed is {F}.")