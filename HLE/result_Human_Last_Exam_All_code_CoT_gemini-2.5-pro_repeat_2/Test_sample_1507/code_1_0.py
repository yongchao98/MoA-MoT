# The problem is to find the minimal number of triangles for triangulating the
# outer-hull of the great icosahedron.

# Step 1: The outer-hull (or convex hull) of the great icosahedron is a
# regular icosahedron, as they share the same 12 vertices.

# Step 2: A regular icosahedron is a Platonic solid, and its surface is
# composed of triangular faces. We need to find the number of these faces.
# By definition, a regular icosahedron has 20 faces.
number_of_faces_on_an_icosahedron = 20

# Step 3: "Triangulating" the surface means dividing it into triangles.
# Since the surface of a regular icosahedron is already made of 20 triangles,
# this is the minimal triangulation possible.
minimal_triangles = number_of_faces_on_an_icosahedron

# Step 4: Print the final answer. The "equation" is simply that the result
# equals the number of faces.
print("The minimal number of triangles required is equal to the number of faces of a regular icosahedron.")
print(f"Result = {minimal_triangles}")