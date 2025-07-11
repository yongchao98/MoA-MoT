# The great icosahedron's outer hull is a regular icosahedron.
# We can use Euler's polyhedron formula (V - E + F = 2) to find the number of faces (F).
# For a regular icosahedron:
# V (Vertices) = 12
# E (Edges) = 30

V = 12
E = 30

# The formula to find the number of faces is F = 2 - V + E
F = 2 - V + E

print("The minimal number of triangles is equal to the number of faces (F) of the regular icosahedron.")
print(f"Using Euler's formula F = 2 - V + E:")
print(f"F = 2 - {V} + {E}")
print(f"The number of triangles needed is: {F}")