# A regular icosahedron has 12 vertices and 30 edges.
# We can use Euler's formula for polyhedra (V - E + F = 2) to find the number of faces (F).
# The formula can be rearranged to F = 2 - V + E.

# Number of vertices in a regular icosahedron
vertices = 12
# Number of edges in a regular icosahedron
edges = 30

# Calculate the number of faces
faces = 2 - vertices + edges

print("The outer hull of a great icosahedron is a regular icosahedron.")
print("The number of triangles for the triangulation is the number of faces of this icosahedron.")
print("Using Euler's formula (F = 2 - V + E):")
print(f"Faces = {2} - {vertices} + {edges}")
print(f"The minimal number of triangles needed is: {faces}")