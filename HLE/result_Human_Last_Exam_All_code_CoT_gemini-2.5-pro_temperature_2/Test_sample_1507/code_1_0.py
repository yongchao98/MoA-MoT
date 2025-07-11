# The great icosahedron is a non-convex polyhedron. Its convex hull (or "outer-hull")
# is a regular icosahedron.
# A regular icosahedron is a Platonic solid made up of 20 identical equilateral triangles.
# Therefore, the minimal number of triangles needed to form its surface is its number of faces.

# Number of faces of a regular icosahedron.
num_triangles = 20

# The problem asks for the minimal number of triangles needed for the triangulation.
# The equation is simply stating this number.
print("The convex hull of the great icosahedron is a regular icosahedron.")
print(f"A regular icosahedron is composed of {num_triangles} triangular faces.")
print(f"Therefore, the minimal number of triangles needed is {num_triangles}.")
