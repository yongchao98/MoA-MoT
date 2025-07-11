# The outer-hull of a great icosahedron is a regular icosahedron.
# A regular icosahedron's surface is already composed of triangles.
# We can find the number of these triangles (faces) using Euler's formula for convex polyhedra: V - E + F = 2.

# For a regular icosahedron:
# V (Vertices) = 12
# E (Edges) = 30
V = 12
E = 30

# We rearrange the formula to solve for F (Faces): F = 2 - V + E
F = 2 - V + E

print("The outer-hull of the great icosahedron is a regular icosahedron.")
print("The number of triangles needed for its triangulation is equal to its number of faces.")
print("Using Euler's formula (V - E + F = 2), we can calculate the number of faces (F).")
print(f"Given V = {V} vertices and E = {E} edges:")
print("The equation to find the number of faces is: F = 2 - V + E")
print(f"Substituting the values: F = 2 - {V} + {E}")
print(f"Result: F = {F}")
print(f"\nThus, the minimal number of triangles is {F}.")
