# The problem is to find the number of triangles on the surface of the
# convex hull of a great icosahedron.

# Step 1: The convex hull of a great icosahedron is a regular icosahedron,
# as they share the same vertex arrangement.

# Step 2: A regular icosahedron is a Platonic solid with 20 faces.
# All of its faces are triangles.

# Step 3: Since the surface is already composed of triangles, the minimal
# number of triangles needed for triangulation is simply its number of faces.
number_of_faces_on_icosahedron = 20

# The final equation is: Number of Triangles = 20
# We will print the number that solves this equation.
print("The minimal number of triangles needed for the triangulation is:")
print(number_of_faces_on_icosahedron)