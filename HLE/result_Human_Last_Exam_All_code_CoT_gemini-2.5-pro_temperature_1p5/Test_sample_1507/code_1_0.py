# Step 1: Define the properties of the regular icosahedron, which is the
# outer-hull of the great icosahedron.

# Number of vertices (V)
V = 12
# Number of triangular faces meeting at each vertex (k)
k = 5
# Number of vertices per triangular face (n_verts_per_face)
n_verts_per_face = 3

# Step 2: Calculate the number of faces (F).
# Each vertex is a part of 'k' faces. If we multiply V * k, we count each
# face 'n_verts_per_face' times (once for each of its vertices).
# So, we divide by 'n_verts_per_face' to get the correct number of faces.
F = (V * k) / n_verts_per_face

# Step 3: Print the explanation and the final equation.
print("The outer-hull of the great icosahedron is a regular icosahedron.")
print("The number of visible triangles is the number of faces on this icosahedron.")
print("\nWe can calculate the number of faces (F) using its vertex properties:")
print(f"Number of Vertices (V) = {V}")
print(f"Faces meeting at each Vertex (k) = {k}")
print(f"Vertices per Face = {n_verts_per_face}")
print("\nThe equation to find the total number of faces is: F = (V * k) / 3")
print("Substituting the values:")
# Using int() for clean output as the result is a whole number.
print(f"F = ({V} * {k}) / {n_verts_per_face} = {int(F)}")
print(f"\nThus, the minimal number of triangles needed is {int(F)}.")