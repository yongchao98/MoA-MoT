# The outer-hull of a great icosahedron is a regular icosahedron.
# We use Euler's formula (V - E + F = 2) to find the number of faces (F).

# A regular icosahedron has 12 vertices.
V = 12
# A regular icosahedron has 30 edges.
E = 30

# Calculate the number of faces (F) using the formula F = 2 - V + E.
F = 2 - V + E

print(f"The outer-hull of the great icosahedron is a regular icosahedron.")
print(f"A regular icosahedron has V = {V} vertices and E = {E} edges.")
print(f"Using Euler's formula V - E + F = 2:")
print(f"{V} - {E} + F = 2")
print(f"Solving for F: F = 2 - {V} + {E}")
print(f"F = {F}")
print(f"\nThe number of visible triangles is {F}.")