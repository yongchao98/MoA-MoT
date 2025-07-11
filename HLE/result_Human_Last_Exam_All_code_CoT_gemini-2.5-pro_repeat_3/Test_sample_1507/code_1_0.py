# The outer-hull of the great icosahedron is a regular icosahedron.
# We can find the number of its faces (F) using Euler's formula: V - E + F = 2
# where V is the number of vertices and E is the number of edges.

# For a regular icosahedron:
# Number of vertices (V)
V = 12
# Number of edges (E)
E = 30

# Using Euler's formula V - E + F = 2, we can solve for F:
# F = 2 - V + E
F = 2 - V + E

print("The outer-hull of the great icosahedron is a regular icosahedron.")
print("To find the number of triangular faces, we use Euler's formula for polyhedra: V - E + F = 2.")
print(f"A regular icosahedron has {V} vertices and {E} edges.")
print("\nPlugging these values into the formula:")
# The prompt asks to output each number in the final equation.
print(f"{V} - {E} + F = 2")
print(f"F = 2 - {V} + {E}")
print(f"F = {F}")

print(f"\nTherefore, the minimal number of triangles needed for triangulating the outer-hull is {F}.")