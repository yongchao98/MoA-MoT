# A regular icosahedron is the convex hull of a great icosahedron.
# We need to find the number of triangular faces on this icosahedron.

# Properties of a regular icosahedron
vertices = 12
edges = 30

# Euler's formula for convex polyhedra is: V - E + F = 2
# We can rearrange it to solve for the number of faces (F): F = 2 - V + E
faces = 2 - vertices + edges

print(f"The convex hull of a great icosahedron is a regular icosahedron.")
print(f"A regular icosahedron has {vertices} vertices (V) and {edges} edges (E).")
print("We use Euler's formula (V - E + F = 2) to find the number of faces (F).")
print("The formula to find F is: F = 2 - V + E")
print(f"So, the number of faces is 2 - {vertices} + {edges} = {faces}")
print(f"\nSince each face of an icosahedron is a triangle, the minimal number of triangles needed is {faces}.")