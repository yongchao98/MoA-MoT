# Step 1: Define the properties of the regular icosahedron (the convex hull).
# It has 12 vertices and 30 edges.
num_vertices = 12
num_edges = 30

# Step 2: Use Euler's formula (V - E + F = 2) to find the number of faces (F).
# Rearranging the formula gives: F = 2 - V + E
num_faces = 2 - num_vertices + num_edges

# Step 3: Print the logic and the final calculation.
# The faces of a regular icosahedron are already triangles, so the number of
# faces is the minimal number of triangles needed for the triangulation.
print("The outer-hull of the great icosahedron is a regular icosahedron.")
print("To find its number of faces (triangles), we use Euler's formula: V - E + F = 2")
print(f"Vertices (V) = {num_vertices}")
print(f"Edges (E) = {num_edges}\n")
print("The equation to find the number of faces is: F = 2 - V + E")
print("Plugging in the numbers:")
# The final part of the request is to show the equation with the numbers.
print(f"2 - {num_vertices} + {num_edges} = {num_faces}")
print(f"\nThe minimal number of triangles needed is {num_faces}.")
