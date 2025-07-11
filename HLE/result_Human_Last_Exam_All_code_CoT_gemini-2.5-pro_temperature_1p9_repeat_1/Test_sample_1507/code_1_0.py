# The convex hull of a great icosahedron is a regular icosahedron.
# We can find its number of faces using its known properties and Euler's formula.

# Number of vertices in a regular icosahedron
vertices = 12
# Number of edges in a regular icosahedron
edges = 30

# Calculate the number of faces using Euler's formula: F = E - V + 2
faces = edges - vertices + 2

print(f"The number of triangles needed is the number of faces of a regular icosahedron.")
print(f"Using Euler's formula (V - E + F = 2), we can find the number of faces (F).")
print(f"Vertices (V) = {vertices}")
print(f"Edges (E) = {edges}")
print(f"The final equation is: F = {edges} - {vertices} + 2")
print(f"The result is: {faces}")