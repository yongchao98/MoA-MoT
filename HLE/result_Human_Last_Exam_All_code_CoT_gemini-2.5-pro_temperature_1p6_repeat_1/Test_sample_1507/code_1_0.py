# The great icosahedron is a non-convex polyhedron. Its outer-hull,
# or convex hull, is a regular icosahedron.
# A regular icosahedron is a Platonic solid with 20 faces.
# Each of these faces is an equilateral triangle.

# The process of "triangulating" a surface means covering it with triangles.
# Since the surface of a regular icosahedron is already made of triangles,
# it is already triangulated in its simplest form.

# Therefore, the minimal number of triangles needed for this triangulation
# is simply the number of faces on a regular icosahedron.

number_of_faces_of_icosahedron = 20

# The final equation is simply the number of faces.
print("The minimal number of triangles needed is the number of faces of a regular icosahedron.")
print(f"The final number is: {number_of_faces_of_icosahedron}")