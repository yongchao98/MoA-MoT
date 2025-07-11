# The great icosahedron is a non-convex polyhedron. Its "outer-hull",
# or the shape visible from the outside, is its convex hull.
# The convex hull of the great icosahedron is a regular icosahedron.

# A regular icosahedron is a platonic solid whose surface is composed
# of triangular faces. The task asks for the minimal number of triangles
# to triangulate this surface. Since the surface is already made of triangles,
# the answer is simply the number of faces of a regular icosahedron.

# The number of faces on a regular icosahedron is a well-known property.
number_of_triangles = 20

# We form the final "equation" by simply stating the result.
print("The outer-hull of the great icosahedron is a regular icosahedron.")
print("The surface of a regular icosahedron consists of 20 triangles.")
print("Final number of triangles needed: {}".format(number_of_triangles))