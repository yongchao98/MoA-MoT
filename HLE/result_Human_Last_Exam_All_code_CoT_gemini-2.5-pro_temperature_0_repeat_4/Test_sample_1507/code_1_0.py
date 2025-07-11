# A great icosahedron's outer-hull (or convex hull) is a regular icosahedron.
# The faces of a regular icosahedron are already triangles.
# We need to find the number of faces of a regular icosahedron.

# We can use Euler's formula for polyhedra: V - E + F = 2
# where V is the number of vertices, E is the number of edges, and F is the number of faces.

# For a regular icosahedron:
# Number of vertices (V)
vertices = 12
# Number of edges (E)
edges = 30

# We can rearrange Euler's formula to solve for F: F = 2 - V + E
faces = 2 - vertices + edges

# Print the equation and the final result
print("The number of visible triangles is the number of faces (F) of the outer-hull (a regular icosahedron).")
print("Using Euler's formula (F = 2 - V + E):")
print(f"F = 2 - {vertices} + {edges}")
print(f"The minimal number of triangles needed is: {faces}")